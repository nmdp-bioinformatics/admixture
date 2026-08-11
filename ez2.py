#!/usr/bin/env python3
"""
Convert Israel subpopulation XLSX (wide *_freq columns) into per-pop CSV files:

  freqs.IL-<POP>.csv

Each output file has columns:
  Haplo,Count,Freq

- Haplo = A~C~B~DRB1
- Freq is normalized to sum to 1 within each population (recommended)
- Count is a pseudo-count = Freq * count_scale (keeps schema consistent with your global files)
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path
from typing import Iterable, List

import pandas as pd


def safe_pop_name(pop: str) -> str:
    """Make population label filename-safe but readable."""
    pop = pop.strip()
    pop = re.sub(r"\s+", "_", pop)
    pop = re.sub(r"[^A-Za-z0-9_\-]+", "", pop)
    return pop


def convert_israel_xlsx_to_freq_files(
    xlsx_path: str,
    out_dir: str = ".",
    keep_pops: Iterable[str] = (
        "ARAB",
        "DRUZE",
        "LIBYA",
        "MAROCCO",
        "TUNISIA",
        "YEMEN",
        "IRAQ",
        "IRAN",
    ),
    prefix: str = "IL-",
    count_scale: float = 1e6,
    sheet_name: str | int | None = None,
    normalize: bool = True,
) -> List[str]:
    """
    Reads the XLSX, selects keep_pops *_freq columns, and writes one file per pop:
      freqs.<prefix><POP>.csv
    Returns list of written file paths.
    """
    out_dir_p = Path(out_dir)
    out_dir_p.mkdir(parents=True, exist_ok=True)

    # Always end up with a DataFrame (not a dict of DataFrames)
    if sheet_name is None:
        # Read first sheet explicitly
        df = pd.read_excel(xlsx_path, sheet_name=0)
    else:
        df = pd.read_excel(xlsx_path, sheet_name=sheet_name)
    
    # Safety: if pandas ever returns a dict anyway, pick the first sheet
    if isinstance(df, dict):
        first_key = next(iter(df.keys()))
        df = df[first_key]


    required_loci = ["A", "C", "B", "DRB1"]
    missing = [c for c in required_loci if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required locus columns in XLSX: {missing}")

    # Build 4-locus haplotype string
    df["Haplo"] = df[required_loci].astype(str).agg("~".join, axis=1)

    # Determine which freq columns exist
    keep_pops = list(keep_pops)
    freq_cols = [f"{p}_freq" for p in keep_pops if f"{p}_freq" in df.columns]
    if not freq_cols:
        raise ValueError(
            "No matching *_freq columns found for keep_pops. "
            f"Requested: {keep_pops}. Available columns include: "
            + ", ".join([c for c in df.columns if c.endswith("_freq")][:20])
            + ("..." if len([c for c in df.columns if c.endswith('_freq')]) > 20 else "")
        )

    # Long format: Haplo, population, Freq
    long = df[["Haplo"] + freq_cols].melt(
        id_vars="Haplo", var_name="population", value_name="Freq"
    )
    long["population"] = long["population"].str.replace("_freq", "", regex=False)
    long["Freq"] = pd.to_numeric(long["Freq"], errors="coerce").fillna(0.0)

    # Keep only positive freqs
    long = long[long["Freq"] > 0].copy()

    written: List[str] = []

    # Write one file per population
    for pop, g in long.groupby("population", sort=True):
        g = g[["Haplo", "Freq"]].copy()

        if normalize:
            s = g["Freq"].sum()
            if s > 0:
                g["Freq"] = g["Freq"] / s

        # Pseudo-count to match schema with your global files
        g["Count"] = g["Freq"] * float(count_scale)

        out = g[["Haplo", "Count", "Freq"]].sort_values("Freq", ascending=False)

        pop_safe = safe_pop_name(pop)
        fname = out_dir_p / f"freqs.{prefix}{pop_safe}.csv"
        out.to_csv(fname, index=False)
        written.append(str(fname))

    return written


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Convert Israel subpopulation XLSX into per-pop freqs.IL-POP.csv files (Haplo,Count,Freq)."
    )
    ap.add_argument("xlsx", help="Input XLSX file (Israel subpopulations wide table)")
    ap.add_argument(
        "--out-dir",
        default=".",
        help="Output directory for freqs.*.csv files (default: current dir)",
    )
    ap.add_argument(
        "--keep",
        default="ARAB,DRUZE,LIBYA,MAROCCO,TUNISIA,YEMEN,IRAQ,IRAN",
        help="Comma-separated population names to export (default: Mediterranean/N. Africa set)",
    )
    ap.add_argument(
        "--prefix",
        default="IL-",
        help="Prefix used in filenames after 'freqs.' (default: IL-)",
    )
    ap.add_argument(
        "--count-scale",
        type=float,
        default=1e6,
        help="Scale factor for pseudo-counts: Count = Freq * count_scale (default: 1e6)",
    )
    ap.add_argument(
        "--sheet",
        default=None,
        help="Optional sheet name or 0-based sheet index (default: first sheet)",
    )
    ap.add_argument(
        "--no-normalize",
        action="store_true",
        help="Do NOT renormalize frequencies within each population after filtering",
    )
    args = ap.parse_args()

    sheet = args.sheet
    if sheet is not None:
        # Allow numeric sheet index passed as string
        try:
            sheet = int(sheet)
        except ValueError:
            pass

    keep_pops = [x.strip() for x in args.keep.split(",") if x.strip()]

    written = convert_israel_xlsx_to_freq_files(
        xlsx_path=args.xlsx,
        out_dir=args.out_dir,
        keep_pops=keep_pops,
        prefix=args.prefix,
        count_scale=args.count_scale,
        sheet_name=sheet,
        normalize=not args.no_normalize,
    )

    for p in written:
        print(f"Wrote: {p}")


if __name__ == "__main__":
    main()

