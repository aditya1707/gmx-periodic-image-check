"""Tests for check_periodic_image.py, using synthetic mindist_pi.xvg files.

These do not require GROMACS: they write small .xvg files by hand and check the
parsing, run-list building, and reported statistics.
"""

import sys
from pathlib import Path

import numpy as np

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


def test_load_min_periodic_skips_headers_and_reads_col2(tmp_path):
    p = tmp_path / "run" / "mindist_pi.xvg"
    write_xvg(p, [1.5, 0.8, 2.0])
    vals = cpi.load_min_periodic(p)
    assert np.allclose(vals, [1.5, 0.8, 2.0])


def test_load_min_periodic_two_column_file(tmp_path):
    p = tmp_path / "run" / "mindist_pi.xvg"
    write_xvg(p, [0.4, 0.6], ncol=2)
    vals = cpi.load_min_periodic(p)
    assert np.allclose(vals, [0.4, 0.6])


def test_analyse_run_statistics(tmp_path):
    # 2 of 5 frames below a 0.3 nm cutoff; closest approach is 0.10 nm.
    p = tmp_path / "run" / "mindist_pi.xvg"
    write_xvg(p, [0.10, 0.25, 0.5, 1.0, 0.9])
    row, series = cpi.analyse_run(tmp_path / "run", "grp", "mindist_pi.xvg", 0.3)
    assert row["n_frames"] == 5
    assert row["closest_approach_nm"] == 0.10
    assert row["frac_within_cutoff"] == 0.4
    assert len(series) == 5


def test_analyse_run_missing_file_returns_none(tmp_path):
    assert cpi.analyse_run(tmp_path / "nope", "grp", "mindist_pi.xvg", 0.3) is None


def test_load_config_inline_comment_and_multiline_runs(tmp_path):
    cfg_file = tmp_path / "config.conf"
    cfg_file.write_text(
        'RUNS="\n'
        '/data/a   sysA\n'
        '/data/b\n'
        '"\n'
        'CUTOFF_NM="0.3"   # inline comment\n'
        'GROUP="Protein"\n'
    )
    cfg = cpi.load_config(cfg_file)
    assert cfg["CUTOFF_NM"] == "0.3"       # inline comment stripped
    assert cfg["GROUP"] == "Protein"
    runs = cpi.build_runs(cfg)
    assert runs == [("/data/a", "sysA"), ("/data/b", "b")]  # label defaults to basename


def test_build_runs_grid_shortcut():
    cfg = {"RUNS": "", "BASE": "/base", "REPLICAS": "1 2", "SYSTEMS": "sysX"}
    runs = cpi.build_runs(cfg)
    assert runs == [("/base/1/sysX", "sysX"), ("/base/2/sysX", "sysX")]
