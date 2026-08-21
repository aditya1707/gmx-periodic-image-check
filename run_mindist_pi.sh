#!/usr/bin/env bash
# =============================================================================
# Generate the periodic-image minimum-distance data for every run in the config.
#
# For each run directory it runs:
#     gmx mindist -pi -f <xtc> -s <tpr> [-n <ndx>] -od <outname>
# which writes a 6-column .xvg: time, min_periodic, max_internal, box_x/y/z.
#
# Usage:
#     ./run_mindist_pi.sh [config_file]      # default: ./config.conf
#     FORCE=1 ./run_mindist_pi.sh            # re-run even where output exists
#
# Existing, non-empty outputs are skipped unless FORCE=1.
# =============================================================================

set -uo pipefail

# --- Locate and load the config ---------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${1:-$SCRIPT_DIR/config.conf}"

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: config file not found: $CONFIG" >&2
    exit 1
fi
# shellcheck source=/dev/null
source "$CONFIG"

# Apply fallbacks so a partially-edited config still runs.
XTC_NAME="${XTC_NAME:-prd.xtc}"
TPR_NAME="${TPR_NAME:-prd.tpr}"
NDX_NAME="${NDX_NAME:-}"
GROUP="${GROUP:-Protein}"
OUTNAME="${OUTNAME:-mindist_pi.xvg}"
FORCE="${FORCE:-0}"
LOG="mindist_pi.log"

# --- gmx must be available ---------------------------------------------------
if ! command -v gmx >/dev/null 2>&1; then
    echo "ERROR: 'gmx' not found on PATH. Load your GROMACS module first." >&2
    exit 1
fi

# --- Build the list of (directory, group label) pairs ------------------------
# Prefer the explicit RUNS list; fall back to the BASE/replica/system grid.
declare -a DIRS GROUPS

runs_stripped="$(printf '%s' "${RUNS:-}" | tr -d '[:space:]')"
if [[ -n "$runs_stripped" ]]; then
    # Explicit list: "directory [group]" per line.
    while read -r dir group _; do
        [[ -z "$dir" || "$dir" == \#* ]] && continue
        DIRS+=("$dir")
        # Default the group label to the directory's own name.
        GROUPS+=("${group:-$(basename "${dir%/}")}")
    done <<< "$RUNS"
else
    # Grid shortcut: BASE/<replica>/<system>, group = system.
    if [[ -n "${BASE:-}" && -n "${REPLICAS:-}" && -n "${SYSTEMS:-}" ]]; then
        for rep in $REPLICAS; do
            for sys in $SYSTEMS; do
                DIRS+=("$BASE/$rep/$sys")
                GROUPS+=("$sys")
            done
        done
    fi
fi

if [[ "${#DIRS[@]}" -eq 0 ]]; then
    echo "ERROR: no runs defined. Fill in RUNS (or the BASE/REPLICAS/SYSTEMS" >&2
    echo "       grid) in $CONFIG." >&2
    exit 1
fi

# --- Process each run --------------------------------------------------------
n_run=0 n_skip=0 n_missing=0 n_fail=0
echo "==== mindist -pi run started $(date) ====" >> "$LOG"

for i in "${!DIRS[@]}"; do
    d="${DIRS[$i]}"
    grp="${GROUPS[$i]}"
    xtc="$d/$XTC_NAME"
    tpr="$d/$TPR_NAME"
    out="$d/$OUTNAME"

    # Skip when a non-empty output already exists (unless forced).
    if [[ -s "$out" && "$FORCE" != "1" ]]; then
        echo "SKIP    $d (output exists; FORCE=1 to redo)"
        n_skip=$((n_skip + 1))
        continue
    fi

    # Required inputs must be present.
    if [[ ! -d "$d" ]]; then
        echo "MISSING $d (directory not found)"
        n_missing=$((n_missing + 1)); continue
    fi
    if [[ ! -f "$xtc" || ! -f "$tpr" ]]; then
        echo "MISSING $d (need $XTC_NAME and $TPR_NAME)"
        n_missing=$((n_missing + 1)); continue
    fi

    # Index file is optional; only pass -n when one is configured and present.
    ndx_arg=()
    if [[ -n "$NDX_NAME" ]]; then
        if [[ -f "$d/$NDX_NAME" ]]; then
            ndx_arg=(-n "$d/$NDX_NAME")
        else
            echo "MISSING $d (NDX_NAME='$NDX_NAME' set but file absent)"
            n_missing=$((n_missing + 1)); continue
        fi
    fi

    echo "RUN     $d  [group: $grp]"
    if echo "$GROUP" | gmx mindist -f "$xtc" -s "$tpr" "${ndx_arg[@]}" \
            -pi -od "$out" >> "$LOG" 2>&1; then
        n_run=$((n_run + 1))
    else
        echo "FAIL    $d (see $LOG)"
        n_fail=$((n_fail + 1))
    fi
done

# --- Summary -----------------------------------------------------------------
echo
echo "Done. ran=$n_run  skipped=$n_skip  missing=$n_missing  failed=$n_fail"
echo "gmx output logged to $LOG"
[[ "$n_fail" -gt 0 ]] && exit 1 || exit 0
