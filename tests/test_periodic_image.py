"""Tests for check_periodic_image.py, using synthetic mindist_pi.xvg files.

These do not require GROMACS: they write small .xvg files by hand and check the
runs-file parsing, the .xvg loading, and the reported statistics.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Make the repo root importable so we can import the script as a module.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import check_periodic_image as cpi  # noqa: E402


def write_xvg(path, min_periodic, ncol=6, box=6.0, max_internal=2.5):
    """Write a synthetic gmx mindist -pi .xvg with the given column-2 values."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("@ title \"minimum periodic distance\"\n")
        f.write("# a comment line\n")
        for i, mp in enumerate(min_periodic):
            if ncol >= 6:
                f.write(f"{i*10.0:.1f} {mp:.4f} {max_internal:.4f} "
                        f"{box:.4f} {box:.4f} {box:.4f}\n")
            else:
                f.write(f"{i*10.0:.1f} {mp:.4f}\n")


def test_read_runs_labels(tmp_path):
    f = tmp_path / "runs.txt"
    f.write_text("# comment\n/data/a   sysA\n/data/b\n\n")
    assert cpi.read_runs(f) == [("/data/a", "sysA"), ("/data/b", "b")]


def test_load_xvg_skips_headers(tmp_path):
    p = tmp_path / "run" / "mindist_pi.xvg"
    write_xvg(p, [1.5, 0.8, 2.0])
    data = cpi.load_xvg(p)
    assert data.shape == (3, 6)
    assert np.allclose(data[:, 1], [1.5, 0.8, 2.0])


def test_analyse_two_column_file(tmp_path):
    p = tmp_path / "run" / "mindist_pi.xvg"
    write_xvg(p, [0.4, 0.6], ncol=2)
    rows, series = cpi.analyse([(str(p.parent), "g")], 0.5, "mindist_pi.xvg")
    assert rows[0]["closest_approach_nm"] == 0.4
    assert np.isnan(rows[0]["max_internal_nm"])   # no cols 3-6 in a 2-col file
    assert len(series) == 1


def test_analyse_statistics_and_informational(tmp_path):
    # 2 of 5 frames below 0.3 nm; closest 0.10 nm; box 4.0, max_internal 3.0.
    p = tmp_path / "run" / "mindist_pi.xvg"
    write_xvg(p, [0.10, 0.25, 0.5, 1.0, 0.9], box=4.0, max_internal=3.0)
    rows, _ = cpi.analyse([(str(p.parent), "g")], 0.3, "mindist_pi.xvg")
    r = rows[0]
    assert r["n_frames"] == 5
    assert r["closest_approach_nm"] == 0.10
    assert r["frac_within_cutoff"] == 0.4
    assert r["max_internal_nm"] == 3.0
    assert r["min_box_nm"] == 4.0
    assert r["min_clearance_nm"] == 1.0        # 4.0 - 3.0


def test_analyse_missing_file(tmp_path, capsys):
    rows, series = cpi.analyse([(str(tmp_path / "nope"), "g")], 0.3, "mindist_pi.xvg")
    assert rows == [] and series == []


def test_pooled_fraction_is_frame_weighted():
    # 100 frames @ 50%  +  900 frames @ 0%  ->  pooled 5%, not the mean 25%.
    df = pd.DataFrame([
        dict(group="s", n_frames=100, frac_within_cutoff=0.50),
        dict(group="s", n_frames=900, frac_within_cutoff=0.00),
    ])
    g = df[df.group == "s"]
    pooled = np.average(g["frac_within_cutoff"], weights=g["n_frames"])
    assert abs(pooled - 0.05) < 1e-9
