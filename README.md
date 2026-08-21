# gmx-periodic-image-check

Check whether a protein — or any chosen group — interacts with its own periodic
image during a GROMACS MD simulation. This is a standard sanity check for a box
that may be too small.

Two steps, driven by a single config file:

1. **`run_mindist_pi.sh`** runs `gmx mindist -pi` on every run you list and
   writes a `mindist_pi.xvg` into each run directory.
2. **`check_periodic_image.py`** reads those files and reports, per run and per
   group, the closest the group ever came to its image and the fraction of
   frames within an interaction cutoff. It writes a CSV, prints a summary, and
   draws a diagnostic time-series plot.

Both read the same **`config.conf`** — the only file you edit.

Built on GROMACS `gmx mindist -pi`; see the
[`gmx mindist` reference](https://manual.gromacs.org/current/onlinehelp/gmx-mindist.html).

## Requirements

- **GROMACS** on your `PATH` (any reasonably recent version; `gmx mindist -pi`
  is long-standing). Needed only for step 1.
- **Python 3.8+** with `numpy` and `pandas`; `matplotlib` is optional (only for
  the plot). See [`requirements.txt`](requirements.txt).

No installation step — clone the repo and run the two scripts.

```bash
git clone <this-repo-url>
cd gmx-periodic-image-check
pip install -r requirements.txt      # numpy, pandas, matplotlib
```

## Use

1. Edit `config.conf`: list your run directories under `RUNS`, set the
   trajectory and topology file names, the `GROUP` to check, and the
   `CUTOFF_NM`.
2. Generate the data (needs GROMACS):

   ```bash
   ./run_mindist_pi.sh
   ```

   Existing, non-empty outputs are skipped; `FORCE=1 ./run_mindist_pi.sh` redoes
   them all. Use a different config with `./run_mindist_pi.sh myconfig.conf`.
3. Summarise:

   ```bash
   python check_periodic_image.py
   ```

   Override the cutoff without editing the config with `--cutoff 1.2`; skip the
   figure with `--no-plot`; point at another config with `--config myconfig.conf`.
   The CSV and plot are written to the directory you run from.

## Layouts

`RUNS` takes one directory per line, so any layout works — a single simulation,
several simulations, or replicas. An optional second word on a line is a group
label that pools related runs (e.g. replicas) in the output. If your data is
laid out as `BASE/<replica>/<system>`, you can instead fill in the
`BASE`/`REPLICAS`/`SYSTEMS` grid shortcut and leave `RUNS` empty.

## Interpreting the result

Both numbers come from **column 2** of the `.xvg` (`min_periodic`) — the minimum
distance between the group and any of its 26 periodic image copies, computed by
GROMACS under full periodic boundary conditions.

- **`closest_approach_nm`** — the closest the group ever came to its image. If
  this stays above your cutoff, the images never interact.
- **`frac_within_cutoff`** — fraction of frames in which that distance fell
  below the cutoff. `0` is the clean result.

The plot is a diagnostic time series: column-2 distance vs frame, with the
cutoff drawn in — useful for seeing *when* (if ever) an image approach happens.

### Why not the max-internal-distance / box check

`gmx mindist -pi` also reports the maximum internal distance and the box
vectors, and a common box-sizing heuristic compares them (a violation when
`box − max_internal < 2·cutoff`). That is a worst-case *a priori* proxy: it
assumes the molecule is a rigid rod aligned with the box vector. For proteins
with large, floppy side chains, `max_internal` is inflated by internal span that
says nothing about image proximity, so the proxy raises false alarms. Column 2
is the direct measurement and is used exclusively here.

## Tests

```bash
pip install pytest
pytest
```

The tests build synthetic `mindist_pi.xvg` files and check the parsing and the
reported statistics; they do not require GROMACS.

## License

MIT — see [`LICENSE`](LICENSE).
