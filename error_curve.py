#!/usr/bin/env python3
"""
Compute NMF reconstruction error for a range of K values and write nmf_error_curve.csv.

Usage:
  python error_curve.py <freqs_dir> [--k-min 2] [--k-max 12]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from nmf import build_matrix, run_nmf


def main():
    ap = argparse.ArgumentParser(description="Compute NMF error curve over a range of K")
    ap.add_argument("freqs_dir", help="Directory with freqs.* files")
    ap.add_argument("--k-min", type=int, default=2, help="Minimum K (default: 2)")
    ap.add_argument("--k-max", type=int, default=12, help="Maximum K (default: 12)")
    ap.add_argument("--min-value", type=float, default=1e-7, help="Min haplotype value filter")
    ap.add_argument("--out", default="nmf_error_curve.csv", help="Output CSV path")
    args = ap.parse_args()

    freqs_dir = Path(args.freqs_dir)
    X, pop_keys, haplo_cols = build_matrix(freqs_dir=freqs_dir, min_value=args.min_value)

    errs = []
    for k in range(args.k_min, args.k_max + 1):
        W, H = run_nmf(X, k=k, random_state=1)
        Xhat = W @ H
        err = np.linalg.norm(X - Xhat, ord="fro")
        errs.append((k, err))
        print(f"  K={k:2d}  error={err:.6f}")

    pd.DataFrame(errs, columns=["K", "error"]).to_csv(args.out, index=False)
    print(f"Wrote: {args.out}")


if __name__ == "__main__":
    main()
