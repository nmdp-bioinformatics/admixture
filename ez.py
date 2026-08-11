#!/usr/bin/env python3
"""
Quick-look script: read an Israel subpopulation XLSX and print summary info.

Usage:
  python ez.py raw/Ezer-Mizion/A~C~B~DRB1.xlsx
"""
import sys
import pandas as pd

if len(sys.argv) < 2:
    print("Usage: python ez.py <xlsx_file>", file=sys.stderr)
    sys.exit(1)

df = pd.read_excel(sys.argv[1])
print(df.head())
print(f"\nShape: {df.shape}")
print(f"Columns: {list(df.columns)}")
