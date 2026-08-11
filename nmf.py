#!/usr/bin/env python3
"""
NMF for haplotype frequency files in your format.

Input: directory containing files like:
  freqs.CY.csv
  freqs.IL-LIBYA.csv
  freqs.US-NAFR.csv  (or freqs.US-NAFR.cav, etc.)

Each file must contain: Haplo,Count,Freq  (extra cols allowed)

We build X (populations x haplotypes) using Freq (default) or Count.
Then run NMF: X ~= W @ H
- W: (n_pops x K) mixture proportions for each original population
- H: (K x n_haplos) component haplotype profiles ("new K populations")

Outputs (in --out-dir):
  W.populations_by_K.csv
  H.K_by_haplotypes.csv
  components.top_haplotypes.csv
  matrix.X.shape.txt
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from sklearn.preprocessing import normalize


def parse_key_from_filename(p: Path) -> str:
    """
    freqs.CY.csv -> CY
    freqs.IL-LIBYA.csv -> IL-LIBYA
    freqs.US-NAFR.cav -> US-NAFR
    """
    name = p.name
    if not name.startswith("freqs."):
        raise ValueError(f"Unexpected file (must start with freqs.): {name}")
    key = name[len("freqs.") :]
    key = re.sub(r"\.[^.]+$", "", key)  # strip last extension
    return key


def read_freq_file(path: Path, use: str) -> pd.DataFrame:
    df = pd.read_csv(path)

    needed = {"Haplo", use}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")

    out = df[["Haplo", use]].copy()
    out[use] = pd.to_numeric(out[use], errors="coerce").fillna(0.0)
    out = out[out[use] > 0].copy()
    return out


def build_matrix(
    freqs_dir: Path,
    use: str = "Freq",
    min_value: float = 0.0,
    renormalize_rows: bool = True,
) -> Tuple[np.ndarray, List[str], List[str]]:
    """
    Returns X (n_pops x n_haplos), pop_keys, haplo_cols
    """
    files = sorted([p for p in freqs_dir.iterdir() if p.is_file() and p.name.startswith("freqs.")])
    if not files:
        raise ValueError(f"No freqs.* files found in {freqs_dir}")

    pop_keys: List[str] = []
    haplo_to_col: Dict[str, int] = {}
    rows_as_dicts: List[Dict[int, float]] = []

    # First pass: gather values and build haplotype index on the fly
    for f in files:
        key = parse_key_from_filename(f)
        df = read_freq_file(f, use=use)

        if min_value > 0:
            df = df[df[use] >= min_value].copy()

        pop_keys.append(key)
        row: Dict[int, float] = {}

        for hap, val in zip(df["Haplo"].astype(str), df[use].astype(float)):
            col = haplo_to_col.get(hap)
            if col is None:
                col = len(haplo_to_col)
                haplo_to_col[hap] = col
            row[col] = row.get(col, 0.0) + float(val)

        rows_as_dicts.append(row)

    n_pops = len(pop_keys)
    n_haplos = len(haplo_to_col)
    X = np.zeros((n_pops, n_haplos), dtype=float)

    for i, row in enumerate(rows_as_dicts):
        for j, v in row.items():
            X[i, j] = v

    # Renormalize each population row to sum to 1 (recommended when using Freq,
    # and still useful when using Count after filtering)
    if renormalize_rows:
        row_sums = X.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        X = X / row_sums

    # Preserve haplotype column order by creation index
    haplo_cols = [""] * n_haplos
    for hap, j in haplo_to_col.items():
        haplo_cols[j] = hap

    return X, pop_keys, haplo_cols


def run_nmf(
    X: np.ndarray,
    k: int,
    random_state: int = 1,
    max_iter: int = 2000,
    alpha_W: float = 0.0,
    alpha_H: float = 0.0,
    l1_ratio: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns W, H with nonnegative entries.
    """
    model = NMF(
        n_components=k,
        init="nndsvdar",
        solver="mu",
        beta_loss="frobenius",
        max_iter=max_iter,
        random_state=random_state,
        alpha_W=alpha_W,
        alpha_H=alpha_H,
        l1_ratio=l1_ratio,
    )
    W = model.fit_transform(X)         # (n_pops x k)
    H = model.components_              # (k x n_haplos)

    # Convert W rows to proportions (sum to 1)
    W = normalize(W, norm="l1", axis=1)

    # Convert H rows to proportions (sum to 1) so each component is a “population profile”
    H = normalize(H, norm="l1", axis=1)

    return W, H


def write_outputs(
    out_dir: Path,
    W: np.ndarray,
    H: np.ndarray,
    pop_keys: List[str],
    haplo_cols: List[str],
    top_n: int = 25,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    k = W.shape[1]
    comp_cols = [f"K{i+1}" for i in range(k)]

    W_df = pd.DataFrame(W, index=pop_keys, columns=comp_cols)
    W_df.index.name = "population"
    W_df.to_csv(out_dir / "W.populations_by_K.csv")

    H_df = pd.DataFrame(H, index=comp_cols, columns=haplo_cols)
    H_df.index.name = "component"
    H_df.to_csv(out_dir / "H.K_by_haplotypes.csv")

    # Top haplotypes per component
    rows = []
    for ci, cname in enumerate(comp_cols):
        weights = H[ci, :]
        top_idx = np.argsort(weights)[::-1][:top_n]
        for rank, j in enumerate(top_idx, start=1):
            rows.append(
                {
                    "component": cname,
                    "rank": rank,
                    "haplo": haplo_cols[j],
                    "weight": float(weights[j]),
                }
            )
    pd.DataFrame(rows).to_csv(out_dir / "components.top_haplotypes.csv", index=False)

    with open(out_dir / "matrix.X.shape.txt", "w", encoding="utf-8") as f:
        f.write(f"X shape: {len(pop_keys)} populations x {len(haplo_cols)} haplotypes\n")


def main():
    ap = argparse.ArgumentParser(description="Run NMF on freqs.*.csv haplotype frequency files")
    ap.add_argument("freqs_dir", help="Directory with freqs.* files")
    ap.add_argument("-k", "--components", type=int, required=True, help="Number of NMF components (K)")
    ap.add_argument("--use", choices=["Freq", "Count"], default="Freq", help="Use Freq or Count as input values")
    ap.add_argument("--min-value", type=float, default=0.0, help="Drop haplotypes below this value before building X")
    ap.add_argument("--no-renorm-rows", action="store_true", help="Do NOT renormalize rows to sum to 1")
    ap.add_argument("--out-dir", default="nmf_out", help="Output directory")
    ap.add_argument("--seed", type=int, default=1, help="Random seed")
    ap.add_argument("--max-iter", type=int, default=2000, help="Max iterations")
    ap.add_argument("--top-n", type=int, default=25, help="Top haplotypes per component to export")

    # Optional regularization knobs (leave at 0 unless you want sparser factors)
    ap.add_argument("--alpha-w", type=float, default=0.0, help="Regularization strength for W")
    ap.add_argument("--alpha-h", type=float, default=0.0, help="Regularization strength for H")
    ap.add_argument("--l1-ratio", type=float, default=0.0, help="0=L2, 1=L1 mix for regularization")

    args = ap.parse_args()

    freqs_dir = Path(args.freqs_dir)
    out_dir = Path(args.out_dir)

    X, pop_keys, haplo_cols = build_matrix(
        freqs_dir=freqs_dir,
        use=args.use,
        min_value=args.min_value,
        renormalize_rows=not args.no_renorm_rows,
    )

    W, H = run_nmf(
        X=X,
        k=args.components,
        random_state=args.seed,
        max_iter=args.max_iter,
        alpha_W=args.alpha_w,
        alpha_H=args.alpha_h,
        l1_ratio=args.l1_ratio,
    )

    write_outputs(out_dir, W, H, pop_keys, haplo_cols, top_n=args.top_n)
    print(f"Wrote outputs to: {out_dir.resolve()}")


if __name__ == "__main__":
    main()

