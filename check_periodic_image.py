#!/usr/bin/env python3
"""Check whether a protein interacts with its own periodic image (GROMACS).

Reads a runs file listing one run directory per line and, for each, uses
``gmx mindist -pi`` output to report the closest the chosen group came to its
periodic image (column 2 of the ``.xvg``) and the fraction of frames within an
interaction cutoff. Column 2 is the verdict; the maximum internal distance and
box vectors (columns 3-6) are reported alongside for information only -- a
worst-case box-sizing view that, for proteins with large floppy side chains, is
inflated by internal span rather than real image proximity.

With ``--run`` it first generates the mindist output by calling gmx; without it,
it summarises existing output files.

Runs file format -- one line per run::

    /path/to/run          optional_group_label

The optional second word pools related runs (e.g. replicas) in the output; if
omitted, the directory name is used. Blank lines and #-comments are ignored.

Note: this is a geometric check on whichever group you select. For a group with
no meaningful internal span (a lone ion or tiny ligand) gmx returns the box
vector rather than a real image distance -- point ``--group`` at the protein.

Usage:
    # generate data (needs gmx), then summarise:
    python check_periodic_image.py runs.txt --run --group Protein --cutoff 0.3
    # summarise existing mindist output, with the diagnostic plot:
    python check_periodic_image.py runs.txt --cutoff 0.3 --plot
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def read_runs(path):
    """Parse the runs file into a list of (directory, group_label) pairs."""
    runs = []
    for ln in Path(path).read_text().splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        parts = ln.split(None, 1)
        directory = parts[0]
        # None marks "no label given" so downstream can fall back appropriately
        # (group -> directory name; plot label -> full path).
        label = parts[1].strip() if len(parts) > 1 else None
        runs.append((directory, label))
    return runs


def generate(runs, args):
    """Run gmx mindist -pi per run, skipping existing non-empty outputs."""
    if shutil.which("gmx") is None:
        sys.exit("ERROR: 'gmx' not found on PATH. Load your GROMACS module first.")

    n_run = n_skip = n_missing = n_fail = 0
    with open("mindist_pi.log", "a") as log:
        for directory, label in runs:
            d = Path(directory)
            out = d / args.outname
            if out.exists() and out.stat().st_size > 0 and not args.force:
                print(f"SKIP    {d} (output exists; --force to redo)")
                n_skip += 1
                continue
            xtc, tpr = d / args.xtc, d / args.tpr
            if not d.is_dir() or not xtc.is_file() or not tpr.is_file():
                print(f"MISSING {d} (need {args.xtc} and {args.tpr})")
                n_missing += 1
                continue
            cmd = ["gmx", "mindist", "-f", str(xtc), "-s", str(tpr),
                   "-pi", "-od", str(out)]
            if args.ndx:
                ndx = d / args.ndx
                if not ndx.is_file():
                    print(f"MISSING {d} (ndx '{args.ndx}' absent)")
                    n_missing += 1
                    continue
                cmd += ["-n", str(ndx)]
            print(f"RUN     {d}  [group: {label or d.name}]")
            log.write(f"\n==== {d} ====\n")
            log.flush()
            r = subprocess.run(cmd, input=args.group + "\n", text=True,
                               stdout=log, stderr=subprocess.STDOUT)
            if r.returncode == 0:
                n_run += 1
            else:
                print(f"FAIL    {d} (see mindist_pi.log)")
                n_fail += 1
    print(f"\nGenerated: ran={n_run} skipped={n_skip} "
          f"missing={n_missing} failed={n_fail}")


def load_xvg(path):
    """Load a gmx mindist -pi .xvg as a 2-D array, skipping @/#/& lines."""
    return np.atleast_2d(np.loadtxt(path, comments=["@", "#", "&"]))


def analyse(runs, cutoff, outname):
    """Compute per-run statistics. Returns (rows, plot_series)."""
    rows, series = [], []
    for directory, label in runs:
        path = Path(directory) / outname
        if not path.exists():
            print(f"  missing: {path}")
            continue
        try:
            data = load_xvg(path)
        except Exception:
            print(f"  unreadable: {path}")
            continue
        if data.size == 0 or data.shape[1] < 2:
            print(f"  unreadable (need >=2 columns): {path}")
            continue

        minper = data[:, 1]
        # Group label for pooling defaults to the directory name; the plot uses
        # the given label, falling back to the full path only when none was set.
        group = label if label else Path(directory).name
        plot_label = label if label else directory
        row = dict(
            group=group,
            run_dir=str(directory),
            n_frames=len(minper),
            closest_approach_nm=round(float(minper.min()), 3),   # the verdict
            frac_within_cutoff=round(float(np.mean(minper < cutoff)), 4),
            max_internal_nm=np.nan,     # informational (columns 3-6 below)
            min_box_nm=np.nan,
            min_clearance_nm=np.nan,
        )
        if data.shape[1] >= 6:
            maxint = data[:, 2]
            box_min = data[:, 3:6].min(axis=1)          # tightest edge per frame
            row["max_internal_nm"] = round(float(maxint.max()), 3)
            row["min_box_nm"] = round(float(box_min.min()), 3)
            # Worst-case clearance between the group's span and a box edge; can
            # go negative when a floppy span exceeds an edge (why column 2 leads).
            row["min_clearance_nm"] = round(float((box_min - maxint).min()), 3)

        rows.append(row)
        series.append((plot_label, minper))
    return rows, series


def report(df, cutoff):
    """Print the per-run table and a per-group verdict + informational block."""
    pd.set_option("display.width", 170)
    pd.set_option("display.max_columns", None)

    print("\n=== Per-run summary ===")
    print(df.to_string(index=False))

    print(f"\n=== Per-group summary (cutoff = {cutoff:g} nm) ===")
    for grp, g in df.groupby("group"):
        closest = g["closest_approach_nm"].min()          # nearest over runs
        # Pool the fraction across runs by frame count, not a plain mean of
        # per-run fractions (runs may differ in length).
        pooled = np.average(g["frac_within_cutoff"], weights=g["n_frames"])
        print(f"\n{grp}  ({len(g)} run(s), {int(g['n_frames'].sum())} frames)")
        print("  verdict -- column 2, actual distance to periodic image:")
        print(f"    closest image approach    : {closest:.3f} nm")
        print(f"    % of frames within {cutoff:g} nm : {pooled * 100:.2f} %")
        if g["max_internal_nm"].notna().any():
            print("  informational -- max internal distance vs box "
                  "(not the verdict):")
            print(f"    max internal distance     : {g['max_internal_nm'].max():.3f} nm")
            print(f"    tightest box edge         : {g['min_box_nm'].min():.3f} nm")
            print(f"    min box-edge clearance    : {g['min_clearance_nm'].min():.3f} nm")


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
    if 1 < len(series) <= 12:
        ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(plot_file, dpi=200)
    print(f"\nPlot written to {plot_file}")


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("runs_file", help="text file listing run directories")
    p.add_argument("--run", action="store_true",
                   help="generate mindist output with gmx before summarising")
    p.add_argument("--force", action="store_true",
                   help="with --run, regenerate even where output exists")
    p.add_argument("--group", default="Protein", help="gmx group to check")
    p.add_argument("--cutoff", type=float, default=1.0,
                   help="interaction cutoff in nm (default 1.0)")
    p.add_argument("--xtc", default="prd.xtc", help="trajectory file name")
    p.add_argument("--tpr", default="prd.tpr", help="run input (.tpr) file name")
    p.add_argument("--ndx", default="", help="index file name (optional)")
    p.add_argument("--outname", default="mindist_pi.xvg",
                   help="mindist output file name, written in each run dir")
    p.add_argument("--csv", default="periodic_image_summary.csv",
                   help="output CSV path")
    p.add_argument("--plot", action="store_true",
                   help="write the diagnostic time-series plot")
    p.add_argument("--plot-file", default="periodic_image_min_distance.png",
                   help="plot output path")
    args = p.parse_args()

    if args.run and not args.tpr.endswith(".tpr"):
        sys.exit(f"ERROR: --tpr '{args.tpr}' is not a .tpr run input file; "
                 f"gmx mindist -pi needs a .tpr.")

    runs = read_runs(args.runs_file)
    if not runs:
        sys.exit(f"No runs found in {args.runs_file}.")

    if args.run:
        generate(runs, args)

    rows, series = analyse(runs, args.cutoff, args.outname)
    if not rows:
        print("\nNo readable output files. Use --run to generate them first.")
        return

    df = pd.DataFrame(rows)
    df.to_csv(args.csv, index=False)
    report(df, args.cutoff)
    print(f"\nCSV written to {args.csv}")

    if args.plot:
        make_plot(series, args.cutoff, args.plot_file)


if __name__ == "__main__":
    main()
