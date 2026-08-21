"""Report whether a protein interacts with its own periodic image.

Reads the ``mindist_pi.xvg`` files produced by ``run_mindist_pi.sh`` (via
``gmx mindist -pi``) for every run listed in ``config.conf`` and reports two
numbers per run and per group:

* closest approach -- the smallest distance the group ever came to its periodic
  image (the minimum of column 2 of the ``.xvg``).
* fraction within cutoff -- the fraction of frames in which that distance fell
  below the interaction cutoff.

Both come from column 2 (``min_periodic``), the minimum distance between the
group and any of its 26 periodic image copies. GROMACS computes this directly
under full periodic boundary conditions, so it is the geometry-exact answer to
"does the group touch its image?". The maximum-internal-distance-vs-box
heuristic is deliberately not used: it is a worst-case box-sizing proxy that is
inflated by large, floppy side chains and does not reflect actual image
proximity.

Outputs a CSV, a short printed summary, and a column-2 time-series diagnostic.

Usage:
    python check_periodic_image.py [--config config.conf] [--cutoff 0.3]
                                   [--csv periodic_image_summary.csv]
                                   [--no-plot] [--plot-file FILE]
"""

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd


# Fallback cutoff (nm) used only if the config has none and --cutoff is unset.
DEFAULT_CUTOFF_NM = 1.0


# --- Config parsing ----------------------------------------------------------
def load_config(path):
    """Parse the shell-style KEY="value" config, including a multi-line RUNS.

    Returns a dict of string values. Handles values quoted on a single line,
    an unquoted single-line value, and a double-quoted value that spans several
    lines (as RUNS does).
    """
    cfg = {}
    lines = Path(path).read_text().splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        i += 1
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, val = line.split("=", 1)
        key = key.strip()
        val = val.strip()
        if val.startswith('"'):
            close = val.find('"', 1)
            if close != -1:
                # Quote closes on the same line; anything after it is a comment.
                cfg[key] = val[1:close]
            else:
                # Open quote with no close: the value continues on later lines.
                buf = [val[1:]]
                while i < len(lines):
                    nxt = lines[i]
                    i += 1
                    if nxt.rstrip().endswith('"'):
                        buf.append(nxt.rstrip()[:-1])
                        break
                    buf.append(nxt)
                cfg[key] = "\n".join(buf)
        else:
            # Unquoted value: drop any trailing inline comment.
            if "#" in val:
                val = val.split("#", 1)[0].strip()
            cfg[key] = val.strip("'")
    return cfg


def build_runs(cfg):
    """Return a list of (directory, group_label) pairs from the config.

    Mirrors run_mindist_pi.sh: the explicit RUNS list wins; otherwise the
    BASE/<replica>/<system> grid is expanded.
    """
    runs = []
    runs_raw = cfg.get("RUNS", "").strip()
    if runs_raw:
        for ln in runs_raw.splitlines():
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            parts = ln.split(None, 1)
            directory = parts[0]
            group = parts[1].strip() if len(parts) > 1 else os.path.basename(
                directory.rstrip("/"))
            runs.append((directory, group))
    else:
        base = cfg.get("BASE", "").strip()
        reps = cfg.get("REPLICAS", "").split()
        syss = cfg.get("SYSTEMS", "").split()
        if base and reps and syss:
            for rep in reps:
                for sys in syss:
                    runs.append((os.path.join(base, rep, sys), sys))
    return runs


# --- Data loading ------------------------------------------------------------
def load_min_periodic(path):
    """Read column 2 (min_periodic, nm) from a gmx mindist -pi .xvg.

    Comment/legend lines (@, #, &) and unparseable lines are skipped. Only the
    first two columns (time, min_periodic) are needed.
    """
    vals = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line[0] in "@#&":
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                vals.append(float(parts[1]))
            except ValueError:
                continue
    return np.asarray(vals, dtype=float)


# --- Per-run analysis --------------------------------------------------------
def analyse_run(directory, group, outname, cutoff):
    """Compute the periodic-image statistics for one run directory.

    Returns (row_dict, min_periodic_series) or None if the file is missing/empty.
    """
    path = Path(directory) / outname
    if not path.exists():
        print(f"  missing: {path}")
        return None
    minper = load_min_periodic(path)
    if minper.size == 0:
        print(f"  unreadable (no data): {path}")
        return None

    row = dict(
        group=group,
        run_dir=str(directory),
        n_frames=len(minper),
        closest_approach_nm=round(float(minper.min()), 3),
        frac_within_cutoff=round(float(np.mean(minper < cutoff)), 4),
    )
    return row, minper


# --- Reporting ---------------------------------------------------------------
def print_report(df, cutoff):
    """Print a per-run table and a per-group closest-approach / fraction summary."""
    pd.set_option("display.width", 160)
    pd.set_option("display.max_columns", None)

    print("\n=== Per-run periodic-image summary ===")
    print(df.to_string(index=False))

    print(f"\n=== Per-group summary (cutoff = {cutoff:g} nm) ===")
    for grp, g in df.groupby("group"):
        closest = g["closest_approach_nm"].min()          # nearest over runs
        mean_frac = g["frac_within_cutoff"].mean()
        print(f"\n{grp}  ({len(g)} run(s))")
        print(f"  closest image approach : {closest:.3f} nm")
        print(f"  % of frames within {cutoff:g} nm   : {mean_frac * 100:.2f} %")


def make_plot(series, cutoff, plot_file):
    """Time-series diagnostic: min-periodic distance (column 2) vs frame."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plot.")
        return

    fig, ax = plt.subplots(figsize=(6, 3.5))
    for label, y in series:
        ax.plot(range(len(y)), y, lw=0.8, label=label)
    ax.axhline(cutoff, ls="--", color="k", lw=1, label=f"cutoff = {cutoff:g} nm")
    ax.set_xlabel("frame")
    ax.set_ylabel("min. distance to periodic image (nm)")
    ax.set_title("Minimum distance to periodic image")
    # A legend helps for a few overlaid runs; skip it when it would be a wall.
    if 1 < len(series) <= 12:
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(plot_file, dpi=200)
    print(f"\nPlot written to {plot_file}")


# --- Entry point -------------------------------------------------------------
def main():
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default=str(script_dir / "config.conf"),
                        help="path to the shared config file")
    parser.add_argument("--cutoff", type=float, default=None,
                        help="interaction cutoff in nm (overrides the config)")
    parser.add_argument("--csv", default="periodic_image_summary.csv",
                        help="output CSV path")
    parser.add_argument("--no-plot", action="store_true", help="skip the plot")
    parser.add_argument("--plot-file", default="periodic_image_min_distance.png",
                        help="output plot path")
    args = parser.parse_args()

    if not Path(args.config).exists():
        parser.error(f"config file not found: {args.config}")
    cfg = load_config(args.config)

    # Cutoff precedence: --cutoff > config > built-in default.
    if args.cutoff is not None:
        cutoff = args.cutoff
    elif cfg.get("CUTOFF_NM", "").strip():
        cutoff = float(cfg["CUTOFF_NM"])
    else:
        cutoff = DEFAULT_CUTOFF_NM

    outname = cfg.get("OUTNAME", "").strip() or "mindist_pi.xvg"
    runs = build_runs(cfg)
    if not runs:
        parser.error("no runs defined in the config (fill in RUNS or the grid).")

    rows, series = [], []
    for directory, group in runs:
        result = analyse_run(directory, group, outname, cutoff)
        if result is not None:
            row, minper = result
            rows.append(row)
            series.append((f"{group}/{Path(directory).name}", minper))

    if not rows:
        print("\nNo readable output files found. Run run_mindist_pi.sh first.")
        return

    df = pd.DataFrame(rows)
    df.to_csv(args.csv, index=False)
    print_report(df, cutoff)
    print(f"\nCSV written to {args.csv}")

    if not args.no_plot:
        make_plot(series, cutoff, args.plot_file)


if __name__ == "__main__":
    main()
