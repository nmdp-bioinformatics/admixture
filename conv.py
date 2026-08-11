#!/usr/bin/env python3
import argparse
import pandas as pd

def collapse_5_to_4(in_csv: str, out_csv: str, renormalize: bool = True) -> None:
    df = pd.read_csv(in_csv)

    # Basic validation
    required = {"Haplo", "Count", "Freq"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    # Keep only required columns; coerce numeric
    df = df[["Haplo", "Count", "Freq"]].copy()
    df["Count"] = pd.to_numeric(df["Count"], errors="coerce").fillna(0.0)
    df["Freq"]  = pd.to_numeric(df["Freq"],  errors="coerce").fillna(0.0)

    # Split and keep first 4 loci
    parts = df["Haplo"].astype(str).str.split("~", expand=True)
    if parts.shape[1] < 4:
        raise ValueError("Haplo column does not contain at least 4 loci delimited by '~'.")

    df["Haplo4"] = parts.iloc[:, 0:4].agg("~".join, axis=1)

    # Sum across DQB1 (and any additional loci)
    out = (
        df.groupby("Haplo4", as_index=False)
          .agg({"Count": "sum", "Freq": "sum"})
          .rename(columns={"Haplo4": "Haplo"})
    )

    # Optional: renormalize frequencies to sum to 1
    if renormalize:
        s = out["Freq"].sum()
        if s > 0:
            out["Freq"] = out["Freq"] / s

    # Nice stable ordering (optional): by descending Freq, then Count
    out = out.sort_values(["Freq", "Count"], ascending=[False, False])

    out.to_csv(out_csv, index=False)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Collapse 5-locus (or longer) haplotypes to first 4 loci and sum across DQB1 variants."
    )
    ap.add_argument("in_csv", help="Input CSV with columns Haplo,Count,Freq (may include extras)")
    ap.add_argument("out_csv", help="Output CSV (4-locus Haplo,Count,Freq)")
    ap.add_argument("--no-renorm", action="store_true",
                    help="Do not renormalize Freq after collapsing")
    args = ap.parse_args()

    collapse_5_to_4(args.in_csv, args.out_csv, renormalize=not args.no_renorm)

