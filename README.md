# gmx-periodic-image-check

Check whether a protein — or any chosen group — interacts with its own periodic
image during a GROMACS MD simulation. A standard sanity check for a box that may
be too small.

One script, `check_periodic_image.py`. It runs `gmx mindist -pi` over a list of
run directories (with `--run`) and reports, per run and per group, the closest
the group ever came to its image and the fraction of frames within an
interaction cutoff. All settings are command-line flags; the only file you write
is a `runs.txt` listing your directories.

Built on GROMACS `gmx mindist -pi`; see the
[`gmx mindist` reference](https://manual.gromacs.org/current/onlinehelp/gmx-mindist.html).

## Requirements

- **GROMACS** on your `PATH` (only for `--run`; `gmx mindist -pi` is
  long-standing).
- **Python 3.8+**. `numpy` and `pandas` are required; `matplotlib` is required
  for `--plot`. Run the script in an environment where these are available;
  `pip install -r requirements.txt` installs them if they are not.

```bash
git clone https://github.com/aditya1707/gmx-periodic-image-check.git
cd gmx-periodic-image-check
pip install -r requirements.txt
```

## Use

1. Write a `runs.txt` — one run directory per line, with an optional group label
   (see [`example_runs.txt`](example_runs.txt)):

   ```
   /data/native/1      native
   /data/native/2      native
   /data/mutant/1      mutant
   /data/one_off_run
   ```

   The optional second word pools related runs (e.g. replicas) in the output; if
   omitted, the directory name is used. Any layout works — one simulation, many,
   with or without replicas.

2. Generate the data and summarise (needs GROMACS):

   ```bash
   python check_periodic_image.py runs.txt --run --group Protein --cutoff 1.0
   ```

   Existing, non-empty outputs are skipped; add `--force` to regenerate.

   To analyse from a given time and **extend** an existing output instead of
   replacing it, pass `-b` (ps). gmx runs from that time into a temporary file
   and only frames later than the existing file's last time are appended, so
   overlap never duplicates frames (a `-b` later than the end leaves a gap):

   ```bash
   python check_periodic_image.py runs.txt --run -b 50000 --group Protein
   ```

   `-e` (ps) sets an end time. Both are passed to `gmx mindist`.

   **Resuming interrupted mindist runs.** If a `gmx mindist` job was killed
   (cluster issue, cancelled), use `--resume`:

   ```bash
   python check_periodic_image.py runs.txt --run --resume --group Protein
   ```

   For each run with an existing output it drops a partial last row (no final
   newline, or fewer columns than the row before), reruns gmx from the last good
   time, and appends the remaining frames. A run that was already complete just
   gets a quick gmx call that adds nothing. Safe to use over a whole runs file.
   It picks its own start time, so it can't be combined with `-b` or `--force`,
   and it requires `--run`. This is about interrupted *mindist* output, not
   unfinished simulations.

3. Re-summarise already-generated output (fast, no GROMACS), with the plot:

   ```bash
   python check_periodic_image.py runs.txt --cutoff 1.0 --plot
   ```

Common flags: `--group` (default `Protein`), `--cutoff` nm (default `1.0`),
`-b`/`--begin` and `-e`/`--end` ps, `--resume`, `--force`, `--xtc`/`--tpr`/`--ndx`/`--outname` file names, `--csv`, `--plot`,
`--plot-file`. The CSV and plot are written to the directory you run from.

## Interpreting the result

The **verdict** comes from **column 2** of the `.xvg` (`min_periodic`) — the
minimum distance between the group and any of its 26 periodic image copies
(26 for `pbc = xyz`; 8 for `pbc = xy`), computed by GROMACS under full periodic
boundary conditions:

- **`closest_approach_nm`** — the closest the group ever came to its image. Above
  the cutoff means the images never interact.
- **`frac_within_cutoff`** — fraction of frames within the cutoff. `0` is clean.

Reported **for information only** (columns 3–6): `max_internal_nm` (the group's
largest internal span), `min_box_nm` (tightest box edge), and `min_clearance_nm`
(`box edge − internal span`, worst case). This is the box-sizing view; it can go
negative for a floppy, extended chain whose span exceeds a box edge without the
protein ever nearing its image — which is exactly why column 2, not this, is the
verdict.

With `--plot` you also get a diagnostic time series: column-2 distance vs frame,
with the cutoff drawn in — useful for seeing *when* (if ever) an approach happens.

## A caveat on the group

This is a geometric check on whichever group you select. For a group with no
meaningful internal span — a lone ion or a tiny ligand — gmx returns the box
vector rather than a real image distance, so the numbers are meaningless. Point
`--group` at the protein (or another extended group).

## Tests

```bash
pip install pytest
pytest
```

The tests build synthetic `mindist_pi.xvg` files and check the parsing and
reported statistics; they do not require GROMACS.

## License

MIT — see [`LICENSE`](LICENSE).
