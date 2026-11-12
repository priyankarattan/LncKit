#!/usr/bin/env bash
# lnckit.sh - launch the LncKit pipeline with all tools on PATH.
#
# LncKit's own code runs on the base Python (3.11, satisfies requires-python
# >=3.10). The bioinformatics tools live in two conda envs and are invoked as
# subprocesses, so they only need to be on PATH -- they do NOT run under this
# interpreter:
#   lnckit         - the main tool set (hisat2, stringtie, cpat, diamond,
#                    samtools 1.18, getorf, cmscan, ...)
#   lnckit_feelnc  - FEELnc only; it pins Perl 5.22 and conflicts with EMBOSS,
#                    so it is isolated in its own env.
#
# lnckit is listed first so its tools win; lnckit_feelnc only fills the FEELnc_*
# gap and its old Perl never shadows anything else.
#
# You run THIS script from the package (by its full path); you pass YOUR working
# directory, and results go THERE - never in the package folder.
#
# Usage (lnckit.sh takes the same arguments as lnckit.py):
#   /path/to/LncKit/lnckit.sh /path/to/my_data          # run the DEFAULT pipeline there
#   /path/to/LncKit/lnckit.sh /path/to/my_data --aligner star --organism plant
#   /path/to/LncKit/lnckit.sh                            # no args -> show the help/options
#
# LncKit is driven entirely by the command line - set only the steps you want to change;
# anything you leave out uses its default. A big run can be launched from a script and
# left on an HPC node (nohup/sbatch/qsub). See every option with:  lnckit.sh --help
#
# The working directory holds your inputs (samples.txt, optionally pipeline.yaml and
# genome files) and receives the lnckit_output/ results folder. The code stays put in
# the package (this script's directory), so one install serves many project folders.
set -euo pipefail

# With NO arguments at all, show the help/options instead of launching a run against
# whatever samples.txt happens to be in the current folder. To run, pass the working
# directory (e.g. `./lnckit.py .`) and/or options.
[ $# -eq 0 ] && set -- --help

# Locate the code dir even when lnckit.sh is reached through a symlink in a project
# folder, so lnckit.py and .lnckit_env are always found next to the real script.
SELF="$(readlink -f "${BASH_SOURCE[0]}" 2>/dev/null || echo "${BASH_SOURCE[0]}")"
HERE="$(cd "$(dirname "$SELF")" && pwd)"

# Find where setup.sh put the environments. Works for conda, mamba OR micromamba, and
# for a shared/read-only HPC base where envs actually land in ~/.conda/envs:
#   1. the ENV_ROOT setup.sh recorded in .lnckit_env, else
#   2. wherever the lnckit env actually exists (checks the usual writable spots).
ENV_ROOT=""
[ -f "$HERE/.lnckit_env" ] && . "$HERE/.lnckit_env" && ENV_ROOT="${LNCKIT_ENV_ROOT:-}"
if [ -z "$ENV_ROOT" ] || [ ! -d "$ENV_ROOT/lnckit" ]; then
    ENV_ROOT=""
    _cands="$HOME/.conda/envs ${MAMBA_ROOT_PREFIX:+$MAMBA_ROOT_PREFIX/envs}"
    command -v conda >/dev/null 2>&1 && _cands="$(conda info --base 2>/dev/null)/envs $_cands"
    _cands="$_cands $HOME/anaconda3/envs $HOME/miniconda3/envs $HOME/miniforge3/envs $HOME/micromamba/envs /opt/conda/envs"
    for c in $_cands; do
        [ -d "$c/lnckit" ] && ENV_ROOT="$c" && break
    done
fi
if [ -z "$ENV_ROOT" ] || [ ! -d "$ENV_ROOT/lnckit" ]; then
    echo "ERROR: LncKit environments not found (looked in: ${ENV_ROOT:-<none>})." >&2
    echo "       Install them first:   bash \"$HERE/setup.sh\"" >&2
    exit 1
fi
LNCKIT_ENV="$ENV_ROOT/lnckit/bin"
FEELNC_ENV="$ENV_ROOT/lnckit_feelnc/bin"
PY="$ENV_ROOT/lnckit/bin/python"   # the lnckit env's own Python runs LncKit itself

# ENV_ROOT is always "<root prefix>/envs"; its parent is the conda/mamba root prefix.
# Export it so any run-time env work (e.g. building the lnckit_py2 env CNCI/PLEK need)
# targets the SAME place setup.sh used - this matters for a micromamba-only server,
# where `micromamba create -n` needs MAMBA_ROOT_PREFIX; harmless for conda.
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$(dirname "$ENV_ROOT")}"

# Only the main `lnckit` env is required. The FEELnc env is optional - if it's not
# installed, just use `classifier: gffcompare` (FEELnc-based classification is skipped).
if [ ! -d "$LNCKIT_ENV" ]; then
    echo "ERROR: main env not found: $LNCKIT_ENV" >&2
    echo "       create it first:  bash \"$HERE/setup.sh\"" >&2
    exit 1
fi
if [ ! -d "$FEELNC_ENV" ]; then
    echo "NOTE: FEELnc env not installed - set 'classifier: gffcompare' in pipeline.yaml." >&2
    FEELNC_ENV=""      # drop it from PATH; FEELnc classification will be skipped
fi

# Tools installed from source/direct downloads (cufflinks, CPC2, PLEK, ...) land
# here; keep it in sync with installer.managed_bin_dir().
MANAGED_BIN="${LNCKIT_HOME:-${XDG_DATA_HOME:-$HOME/.local/share}}/lnckit/bin"

export PATH="$LNCKIT_ENV${FEELNC_ENV:+:$FEELNC_ENV}:$MANAGED_BIN:$PATH"

# Make the lnckit package importable without needing `pip install` - the code is
# right here next to this script.
export PYTHONPATH="$HERE:${PYTHONPATH:-}"

# Shared reference databases (UniProt, Rfam) live here by default - inside the
# package's databases/ folder - so they are downloaded ONCE and reused by every
# project. Override by exporting LNCKIT_DB_DIR to your own location.
export LNCKIT_DB_DIR="${LNCKIT_DB_DIR:-$HERE/databases}"

# Parse options:
#   -o/--output <dir>   results directory
#   --from  <step>      START the pipeline at this step (skip everything before it)
#   --until <step>      STOP after this step
# Steps: fetch merge qc trim index align assemble compare filter_length filter_exon
#        filter_expr filter_orf coding protein housekeeping identify classify report
# e.g. you already have SAM/BAM ->  lnckit.sh <dir> --from align
# Remaining args: an optional WORKING directory (where inputs live) + a sample sheet.
OUTPUT_DIR=""
POS=()
while [ $# -gt 0 ]; do
    case "$1" in
        -o|--output)          OUTPUT_DIR="$2"; shift 2;;
        --output=*)           OUTPUT_DIR="${1#*=}"; shift;;
        -c|--config)          export LNCKIT_CONFIG="$2"; shift 2;;
        --config=*)           export LNCKIT_CONFIG="${1#*=}"; shift;;
        --organism|-org)      export LNCKIT_ORGANISM="$2"; shift 2;;
        --organism=*)         export LNCKIT_ORGANISM="${1#*=}"; shift;;
        --strandness|--strand) export LNCKIT_STRANDNESS="$2"; shift 2;;   # RF | FR | none
        --strandness=*)       export LNCKIT_STRANDNESS="${1#*=}"; shift;;
        --from|--from-step)   export LNCKIT_FROM_STEP="$2"; shift 2;;
        --from=*)             export LNCKIT_FROM_STEP="${1#*=}"; shift;;
        --until|--until-step) export LNCKIT_UNTIL_STEP="$2"; shift 2;;
        --until=*)            export LNCKIT_UNTIL_STEP="${1#*=}"; shift;;
        *)                    POS+=("$1"); shift;;
    esac
done
set -- "${POS[@]+"${POS[@]}"}"

# Working directory (where inputs live and, by default, results are written):
#   - first positional arg, if it is an existing directory -> use it (and consume it)
#   - otherwise the current directory
WORKDIR="$PWD"
if [ "${1:-}" ] && [ -d "$1" ]; then
    WORKDIR="$(cd "$1" && pwd)"
    shift
fi
# Resolve --config to an ABSOLUTE path while we're still in the invocation dir (before the
# cd below), so a relative path like 'my_config.yaml' still works after cd. If it
# isn't found relative to here, try next to the package (HERE) - that's where presets ship.
if [ -n "${LNCKIT_CONFIG:-}" ]; then
    if [ -f "$LNCKIT_CONFIG" ]; then
        LNCKIT_CONFIG="$(cd "$(dirname "$LNCKIT_CONFIG")" && pwd)/$(basename "$LNCKIT_CONFIG")"
    elif [ -f "$HERE/$LNCKIT_CONFIG" ]; then
        LNCKIT_CONFIG="$HERE/$LNCKIT_CONFIG"
    else
        echo "ERROR: --config file not found: $LNCKIT_CONFIG" >&2; exit 2
    fi
    export LNCKIT_CONFIG
    echo "Config          : $LNCKIT_CONFIG"
fi
cd "$WORKDIR"

# If -o/--output was given, results go THERE (resolved to an absolute path) instead
# of WORKDIR/lnckit_output. Passed to lnckit.py via LNCKIT_OUTPUT_DIR.
RESULTS="$WORKDIR/lnckit_output"
if [ -n "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
    OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
    export LNCKIT_OUTPUT_DIR="$OUTPUT_DIR"
    RESULTS="$OUTPUT_DIR"
fi

# Stay quiet for a plain `--help`/`-h` so the option list isn't buried under setup lines.
_SHOW_HELP=0
for _a in "$@"; do case "$_a" in -h|--help) _SHOW_HELP=1;; esac; done
if [ "$_SHOW_HELP" = 0 ]; then
    echo "LncKit package : $HERE"
    echo "Working dir     : $WORKDIR"
    echo "Results go to   : $RESULTS"
fi
# Tell lnckit.py the environment is already set up (tools on PATH), so that when
# lnckit.py is invoked directly it doesn't bounce back here in a loop.
export _LNCKIT_BOOTSTRAPPED=1
exec "$PY" "$HERE/lnckit.py" "$@"
