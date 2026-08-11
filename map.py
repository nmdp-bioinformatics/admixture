#!/usr/bin/env python3
"""
Map Mediterranean / North African populations from freqs.*.csv files.

Input files (examples):
  freqs.CY.csv
  freqs.IL-LIBYA.csv
  freqs.US-NAFR.cav  (treated as US-NAFR; plotted over North Africa as "unspecified")

Each CSV must contain a numeric 'Count' column.

Output:
  med_pop_map.png
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


def parse_population_key(path: Path) -> str:
    """
    Extract key from filenames like:
      freqs.CY.csv -> CY
      freqs.IL-LIBYA.csv -> IL-LIBYA
      freqs.US-NAFR.cav -> US-NAFR
    """
    name = path.name
    if not name.startswith("freqs."):
        raise ValueError(f"Unexpected filename (does not start with freqs.): {name}")
    key = name[len("freqs.") :]
    # strip final extension (.csv, .cav, etc.)
    key = re.sub(r"\.[^.]+$", "", key)
    return key


def sum_count(csv_path: Path) -> float:
    df = pd.read_csv(csv_path)
    if "Count" not in df.columns:
        raise ValueError(f"Missing Count column in {csv_path}")
    c = pd.to_numeric(df["Count"], errors="coerce").fillna(0.0).sum()
    return float(c)


def load_naturalearth_lowres():
    """
    Loads Natural Earth via GeoPandas if available.
    """
    try:
        import geopandas as gpd
    except ImportError as e:
        raise SystemExit(
            "This script needs geopandas. Install with:\n"
            "  pip install geopandas\n"
            "On macOS you may also need:\n"
            "  pip install pyproj shapely fiona"
        ) from e

    world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
    return world


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("freqs_dir", help="Folder containing freqs.* files")
    ap.add_argument("--out", default="med_pop_map.png", help="Output PNG filename")
    ap.add_argument("--bbox", default="-12,45,25,48",
                    help="Map bounding box lon_min,lon_max,lat_min,lat_max (default: Mediterranean)")
    ap.add_argument("--min-size", type=float, default=40.0, help="Minimum bubble size")
    ap.add_argument("--max-size", type=float, default=1200.0, help="Maximum bubble size")
    args = ap.parse_args()

    lon_min, lon_max, lat_min, lat_max = [float(x) for x in args.bbox.split(",")]

    freqs_dir = Path(args.freqs_dir)
    files = sorted([p for p in freqs_dir.iterdir() if p.is_file() and p.name.startswith("freqs.")])

    if not files:
        raise SystemExit(f"No files starting with 'freqs.' found in {freqs_dir}")

    rows = []
    for f in files:
        key = parse_population_key(f)
        total = sum_count(f)
        rows.append({"key": key, "count_sum": total})

    pops = pd.DataFrame(rows)

    # Derive "country code" from key (prefix before first '-'; else whole key)
    pops["country"] = pops["key"].str.split("-", n=1).str[0]

    # Special handling: US-NAFR is "unspecified North African" — place over N. Africa
    # You can move this coordinate if you want (currently near Tunisia/Algeria region).
    special_points = {
        "US-NAFR": {"lon": 9.0, "lat": 34.0},  # “North Africa” anchor
    }

    # Load map
    world = load_naturalearth_lowres()

    # Crop to bbox
    world_crop = world.cx[lon_min:lon_max, lat_min:lat_max].copy()

    # Compute centroids for countries in bbox (GeoPandas geometry)
    # We'll use ISO_A2 codes are not reliably present in this dataset; use 'iso_a3' instead.
    # So map your 2-letter codes -> iso_a3 via a small lookup.
    iso2_to_iso3 = {
        "CY": "CYP",
        "HR": "HRV",
        "ES": "ESP",
        "FR": "FRA",
        "GR": "GRC",
        "IT": "ITA",
        "TR": "TUR",
        "IL": "ISR",
        "US": "USA",
    }

    pops["iso3"] = pops["country"].map(iso2_to_iso3)

    # Find centroid lon/lat for each iso3 from Natural Earth
    centroids = (
        world[["iso_a3", "geometry"]]
        .dropna()
        .drop_duplicates(subset=["iso_a3"])
        .copy()
    )
    centroids["centroid"] = centroids.geometry.centroid
    centroids["lon"] = centroids["centroid"].x
    centroids["lat"] = centroids["centroid"].y

    pops = pops.merge(centroids[["iso_a3", "lon", "lat"]], left_on="iso3", right_on="iso_a3", how="left")

    # Override with special points
    for k, pt in special_points.items():
        pops.loc[pops["key"] == k, ["lon", "lat"]] = (pt["lon"], pt["lat"])

    # Drop anything we couldn't locate
    missing_loc = pops[pops["lon"].isna() | pops["lat"].isna()]
    if len(missing_loc) > 0:
        print("Warning: Could not locate these keys on map (missing iso mapping):")
        for k in missing_loc["key"].tolist():
            print("  -", k)
        pops = pops.dropna(subset=["lon", "lat"])

    # Scale bubble sizes by count_sum (log-ish scaling to avoid one huge point)
    counts = pops["count_sum"].values
    if counts.max() <= 0:
        sizes = [args.min_size] * len(pops)
    else:
        # sqrt scaling tends to look good
        scaled = (counts / counts.max()) ** 0.5
        sizes = args.min_size + scaled * (args.max_size - args.min_size)

    # Plot
    fig, ax = plt.subplots(figsize=(11, 7))

    # base map
    world_crop.plot(ax=ax, linewidth=0.6, edgecolor="black", color="white")

    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title("Phoenician Route: Mediterranean / North African Populations\n(bubble size ~ summed Count)")

    ax.scatter(pops["lon"], pops["lat"], s=sizes, alpha=0.55)

    # labels
    def fmt_count(x: float) -> str:
        # readable compact formatting
        if x >= 1e6:
            return f"{x/1e6:.2f}M"
        if x >= 1e3:
            return f"{x/1e3:.1f}k"
        return f"{x:.0f}"

    for _, r in pops.iterrows():
        label = f"{r['key']} ({fmt_count(r['count_sum'])})"
        ax.text(r["lon"], r["lat"], label, fontsize=8, ha="left", va="bottom")

    # Footnote for special points
    if any(pops["key"].isin(special_points.keys())):
        ax.text(
            0.01, 0.01,
            "Note: US-NAFR plotted as an 'unspecified North Africa' anchor point.",
            transform=ax.transAxes,
            fontsize=8,
            ha="left", va="bottom"
        )

    plt.tight_layout()
    plt.savefig(args.out, dpi=200)
    print(f"Wrote: {args.out}")


if __name__ == "__main__":
    main()

