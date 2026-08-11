# Admixture via NMF

Non-negative Matrix Factorization (NMF) applied to multi-locus HLA haplotype
frequency data to decompose populations into admixture components.

## Overview

Given a set of population-level haplotype frequency files, this pipeline:

1. Builds a matrix **X** (populations × haplotypes)
2. Factorizes **X ≈ W × H** using NMF for a range of K values
3. Visualizes results as stacked-bar admixture plots, pie-chart maps, choropleth
   maps for individual haplotypes, and elbow/error curves

## Directory layout

```
.
├── raw/              # (gitignored) raw input data files
├── freqs/            # (gitignored) processed per-population CSV files
├── nmf_out_k*/       # (gitignored) NMF output for each K
├── nmf.py            # core NMF decomposition
├── error_curve.py    # compute reconstruction error vs K
├── run.sh            # run NMF for K=2..12
├── map.py            # population bubble map
├── ez2.py            # convert Israel-subpop XLSX → per-pop CSV
├── ez.py             # quick-look at an XLSX
├── conv.py           # collapse 5-locus haplotypes to 4 loci
├── map.R             # pie-chart map + stacked bars for all K
├── nmfplot.R         # stacked-bar plot for a single K
├── elbow.R           # elbow/error-reduction plot
├── chloropleth.R     # choropleth map for a single haplotype
├── gp.txt            # ISO-2 country codes used in global data
└── requirements.txt  # Python dependencies
```

## Setup

### Python

```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

### R

Install the following R packages (CRAN):

```r
install.packages(c(
  "readr", "dplyr", "tidyr", "ggplot2", "stringr", "forcats",
  "sf", "rnaturalearth", "scatterpie", "scales",
  "readxl", "png", "grid", "gridExtra"
))
```

## Preparing data

Input frequency files must be CSVs with at minimum: `Haplo`, `Count`, `Freq`.
Place them in the `freqs/` directory with filenames like:

```
freqs.<POPULATION_KEY>.csv
```

For example: `freqs.CY.csv`, `freqs.IL-LIBYA.csv`, `freqs.US-NAFR.csv`

### Converting source data

If your source data is in a different format, use the conversion scripts:

```bash
# Israel subpopulation XLSX → per-pop CSVs
python ez2.py raw/<source>.xlsx --out-dir freqs/

# Collapse 5-locus to 4-locus haplotypes
python conv.py <input.csv> <output.csv>
```

## Running the NMF pipeline

```bash
# Run NMF for K=2 through K=12 and generate error curve
./run.sh freqs/

# Or run for a single K value:
python nmf.py freqs/ -k 6 --out-dir nmf_out_k6
```

## Visualization

```bash
# Stacked-bar admixture plot for a given K
Rscript nmfplot.R 6

# Pie-chart maps + stacked bars for all K values (K=2..12)
Rscript map.R

# Elbow plot (reconstruction error vs K)
Rscript elbow.R 12

# Choropleth for a specific haplotype
Rscript chloropleth.R freqs/ "A*33:01~C*08:02~B*14:02~DRB1*01:02"

# Population bubble map
python map.py freqs/
```

## Output

All generated figures (PNG/PDF) and NMF output directories are gitignored.
They will be created in the working directory when you run the scripts.
