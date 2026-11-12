#!/usr/bin/env bash
# setup.sh - install everything LncKit needs, on ANY Linux system (server, workstation,
# or HPC). No pre-installed conda required: if none is found this downloads micromamba
# (a single self-contained binary, no root needed) and uses that.
#
# It creates three isolated environments:
#   lnckit           - LncKit's own Python + all the bio tools (hisat2, stringtie,
#                      cpat, cpc2, diamond, infernal, getorf, gffcompare, ...)
#   lnckit_feelnc    - FEELnc (old Perl 5.22) + the Berkeley-DB compat fix it needs
#   lnckit_lncfinder - LncFinder deps (R) + LncFinder from CRAN
#
# CNCI and PLEK (not conda packages) are built from source on the first ./lnckit.sh.
# Idempotent: safe to re-run. Usage:  bash setup.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---------------------------------------------------------------------------
# 0. Find a conda-family package manager, or bootstrap micromamba (no root).
#    Sets: PKG (create command prefix) and ENV_ROOT (where envs live).
# ---------------------------------------------------------------------------
# conda/mamba put NEW envs in the first WRITABLE envs dir: the base's envs/ if you can
# write there, else ~/.conda/envs - the usual case on a shared/read-only HPC install
# (e.g. /cm/shared/apps/anaconda3). Detect that so we look where envs REALLY land.
_conda_envs_dir() {
    # Use <base>/envs if the base is writable (create the dir if it doesn't exist yet),
    # else ~/.conda/envs. Testing base/envs directly fails when it doesn't exist yet
    # (fresh conda), which wrongly sent envs to ~/.conda even when the base was writable.
    local base="${1:-}"
    if [ -n "$base" ] && mkdir -p "$base/envs" 2>/dev/null && [ -w "$base/envs" ]; then
        echo "$base/envs"
    else
        mkdir -p "$HOME/.conda/envs" 2>/dev/null || true
        echo "$HOME/.conda/envs"
    fi
}
resolve_pkg() {
    # FAST BY DEFAULT - no flags, no extra commands needed. LncKit resolves its large
    # bioconda environment with a libmamba-based solver (micromamba/mamba), so a plain
    # `bash setup.sh` finishes in minutes instead of the hours the old classic conda solver
    # can take - the usual cause of a setup that "runs all night and never loads". If no fast
    # manager is installed, the bootstrap block below provisions micromamba automatically:
    # from conda-forge when a conda is already present (reliable on HPCs whose outbound
    # access is limited to the conda channels), otherwise by direct download. Nothing for the
    # user to install. Envs live under the manager's envs dir; lnckit.sh finds them there.
    #   LNCKIT_USE_CONDA=1 forces the system conda/mamba instead (may be slow on an old solver).
    if [ "${LNCKIT_USE_CONDA:-0}" = "1" ]; then
        if command -v conda >/dev/null 2>&1; then
            PKG="conda"; ENV_ROOT="$(_conda_envs_dir "$(conda info --base 2>/dev/null)")"
        elif command -v mamba >/dev/null 2>&1; then
            PKG="mamba"; ENV_ROOT="$(_conda_envs_dir "")"
        else
            return 1
        fi
        mkdir -p "$ENV_ROOT" 2>/dev/null || true
        return 0
    fi
    # Default: use an already-installed FAST manager if present (micromamba, then mamba).
    if command -v micromamba >/dev/null 2>&1; then
        export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
        PKG="micromamba"; ENV_ROOT="$MAMBA_ROOT_PREFIX/envs"
        mkdir -p "$ENV_ROOT" 2>/dev/null || true
        return 0
    fi
    if command -v mamba >/dev/null 2>&1; then
        PKG="mamba"; ENV_ROOT="$(_conda_envs_dir "")"
        mkdir -p "$ENV_ROOT" 2>/dev/null || true
        return 0
    fi
    return 1    # no fast manager installed -> bootstrap micromamba (below)
}

if ! resolve_pkg; then
    echo "[bootstrap] no fast package manager found - provisioning micromamba (no root)..."
    BIN="$HOME/.local/bin"; mkdir -p "$BIN"
    MM=""
    # Pick the micromamba build for THIS OS/arch so the download works beyond linux-64 (macOS
    # and ARM/PowerPC Linux too). micro.mamba.pm serves these platform keys:
    case "$(uname -s)-$(uname -m)" in
        Linux-x86_64)              MM_PLAT="linux-64" ;;
        Linux-aarch64|Linux-arm64) MM_PLAT="linux-aarch64" ;;
        Linux-ppc64le)             MM_PLAT="linux-ppc64le" ;;
        Darwin-x86_64)             MM_PLAT="osx-64" ;;
        Darwin-arm64)              MM_PLAT="osx-arm64" ;;
        *)                         MM_PLAT="linux-64" ;;   # best-effort default
    esac
    # 1) Direct static-binary download (fast when micro.mamba.pm is reachable). The final
    #    `--version` run-check rejects a wrong-arch/corrupt binary so we fall through to the
    #    conda-forge fallback instead of accepting something that cannot actually execute.
    if curl -fsSL "https://micro.mamba.pm/api/micromamba/$MM_PLAT/latest" -o "$BIN/_mm.tar.bz2" 2>/dev/null \
       && tar -xjf "$BIN/_mm.tar.bz2" -C "$BIN" --strip-components=1 bin/micromamba 2>/dev/null \
       && [ -x "$BIN/micromamba" ] && "$BIN/micromamba" --version >/dev/null 2>&1; then
        MM="$BIN/micromamba"; echo "[bootstrap] micromamba downloaded ($MM_PLAT)"
    fi
    rm -f "$BIN/_mm.tar.bz2"
    # 2) Fallback: if a conda/mamba exists, install micromamba FROM conda-forge. This works on
    #    HPCs where micro.mamba.pm is blocked but the conda channels are allow-listed. Solving
    #    one small package is quick even with the classic solver.
    if [ -z "$MM" ] && { command -v conda >/dev/null 2>&1 || command -v mamba >/dev/null 2>&1; }; then
        echo "[bootstrap] direct download unavailable - installing micromamba from conda-forge..."
        MMPFX="$HOME/.lnckit-micromamba"
        if command -v mamba >/dev/null 2>&1; then
            mamba create -y -p "$MMPFX" -c conda-forge micromamba >/dev/null 2>&1 || true
        else
            conda create -y -p "$MMPFX" -c conda-forge micromamba >/dev/null 2>&1 || true
        fi
        if [ -x "$MMPFX/bin/micromamba" ]; then
            cp "$MMPFX/bin/micromamba" "$BIN/micromamba" && MM="$BIN/micromamba" \
                && echo "[bootstrap] micromamba installed from conda-forge"
        fi
    fi
    if [ -n "$MM" ]; then
        export PATH="$BIN:$PATH"
        export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
        micromamba shell init -s bash -r "$MAMBA_ROOT_PREFIX" >/dev/null 2>&1 || true
        resolve_pkg
    elif command -v conda >/dev/null 2>&1 || command -v mamba >/dev/null 2>&1; then
        echo "[bootstrap] could not get micromamba - falling back to system conda/mamba (slower)." >&2
        LNCKIT_USE_CONDA=1 resolve_pkg || { echo "ERROR: no usable package manager found." >&2; exit 1; }
    else
        echo "ERROR: no conda/mamba present and micromamba could not be downloaded." >&2
        echo "       Re-run on a node with internet access, or install conda/micromamba first." >&2
        exit 1
    fi
fi
echo "[env] package manager: $PKG    env root: $ENV_ROOT"

cat <<'EOF'

============================================================
  LncKit setup - downloading & installing the tools.
  Time on a fresh machine: roughly 20-40 min (mostly
  downloads; the R/LncFinder env is the slowest step).
  Tip: `LNCKIT_SKIP_LNCFINDER=1 bash setup.sh` skips that
       env and roughly halves the time.
  (The first ./lnckit.sh then builds CNCI/PLEK and downloads
   the UniProt/Rfam databases - another ~10-30 min, once.)
============================================================
EOF

# unified helpers ------------------------------------------------------------
# conda_env_path <name>: the REAL path of a named env, straight from conda (empty if
# absent). `conda env list` is authoritative - it knows the actual envs_dir a named env
# landed in, which can differ from our guess on a read-only shared/HPC base.
conda_env_path() { conda env list 2>/dev/null | awk -v n="$1" '$1==n{print $NF; exit}'; }
# A named env may live in a DIFFERENT envs_dir than we predicted (shared base vs
# ~/.conda), so for conda/mamba ask conda authoritatively; micromamba uses ENV_ROOT.
env_exists() {
    if [ "$PKG" != "micromamba" ] && command -v conda >/dev/null 2>&1; then
        conda env list 2>/dev/null | awk -v n="$1" '$1==n{f=1} END{exit !f}'
    else
        [ -d "$ENV_ROOT/$1" ]
    fi
}
# Prefer `conda env create` (stable exit codes) whenever conda exists - `mamba env
# create` can return non-zero even on success, which trips `set -e`. micromamba (no
# conda) uses its own create.
create_env() {   # $1 = env yml file
    if [ "$PKG" = "micromamba" ]; then
        micromamba create -y -f "$1"
    elif command -v conda >/dev/null 2>&1; then
        conda env create -f "$1"
    else
        mamba env create -f "$1"
    fi
}
create_prefix() {  # $@ = args after `create -y -p ...`
    if [ "$PKG" = "micromamba" ]; then
        micromamba create -y "$@"
    elif command -v conda >/dev/null 2>&1; then
        conda create -y "$@"
    else
        mamba create -y "$@"
    fi
}
# FEELnc solves faster/more reliably via the channel-qualified spec `bioconda::feelnc`
# than via the env yml (the qualifier constrains the search space) - per the FEELnc docs.
create_feelnc_env() {
    if [ "$PKG" = "micromamba" ]; then
        micromamba create -y -n lnckit_feelnc -c conda-forge -c bioconda feelnc
    elif command -v conda >/dev/null 2>&1; then
        conda create -y -n lnckit_feelnc bioconda::feelnc
    else
        mamba create -y -n lnckit_feelnc bioconda::feelnc
    fi
}

# ---- 1. main tool env ------------------------------------------------------
if env_exists lnckit; then echo "[skip]  env 'lnckit' already exists"
else echo "[create] env 'lnckit'"; create_env "$HERE/envs/lnckit.yml"; fi

# RE-POINT ENV_ROOT at where 'lnckit' ACTUALLY landed. Our _conda_envs_dir() only PREDICTS
# the envs_dir; conda may place a named env elsewhere (a read-only shared/HPC base falls
# back to ~/.conda, etc.). If we guessed wrong, every $ENV_ROOT/lnckit/bin/python below
# would be "no such file or directory" (the tool pre-install / MultiQC / CPAT steps). So
# ask conda for the real path and correct ENV_ROOT before anything uses it.
if [ "$PKG" != "micromamba" ] && command -v conda >/dev/null 2>&1; then
    _real="$(conda_env_path lnckit)"
    if [ -n "$_real" ] && [ -d "$_real" ] && [ "$(dirname "$_real")" != "$ENV_ROOT" ]; then
        echo "[env]   'lnckit' is at $_real - re-pointing env root (was $ENV_ROOT)"
        ENV_ROOT="$(dirname "$_real")"
    fi
fi

# ---- 2. FEELnc env ---------------------------------------------------------
# FEELnc's old-Perl/BioPerl stack is the slowest to solve. Skip it with
# LNCKIT_SKIP_FEELNC=1 and use the gffcompare classifier instead.
if [ "${LNCKIT_SKIP_FEELNC:-0}" = "1" ]; then
    echo "[skip]  FEELnc env (LNCKIT_SKIP_FEELNC=1) - use 'classifier: gffcompare'"
elif env_exists lnckit_feelnc; then echo "[skip]  env 'lnckit_feelnc' already exists"
else echo "[create] env 'lnckit_feelnc'  (conda create -n lnckit_feelnc bioconda::feelnc)"
     # FEELnc is OPTIONAL (its old Perl/BioPerl stack is the hardest to solve). A build
     # failure must NOT abort setup - fall back to the gffcompare classifier / builtin filter.
     create_feelnc_env || echo "[warn]  FEELnc env failed to build - skipping it; set 'classifier: gffcompare' and 'candidate_filter: builtin'"
fi

# ---- 3. LncFinder env (deps from conda, LncFinder from CRAN) ---------------
if [ "${LNCKIT_SKIP_LNCFINDER:-0}" = "1" ]; then
    echo "[skip]  LncFinder env (LNCKIT_SKIP_LNCFINDER=1)"
elif env_exists lnckit_lncfinder; then
    echo "[skip]  env 'lnckit_lncfinder' already exists"
else
    echo "[create] env 'lnckit_lncfinder' (deps)"
    # LncFinder is OPTIONAL (one coding-ensemble backend). A deps-env build failure must NOT
    # abort setup - the ensemble simply drops LncFinder from the vote.
    create_env "$HERE/envs/lnckit_lncfinder.yml" || echo "[warn]  LncFinder deps env failed to build - the ensemble will skip LncFinder"
fi
LF_RSCRIPT="$ENV_ROOT/lnckit_lncfinder/bin/Rscript"
if [ -x "$LF_RSCRIPT" ] && ! "$LF_RSCRIPT" -e 'quit(status=!requireNamespace("LncFinder",quietly=TRUE))' 2>/dev/null; then
    echo "[fix]   installing LncFinder (prebuilt binary when available; parallel compile otherwise)"
    # Speed: try Posit Public Package Manager BINARY packages for this Linux distro first
    # (downloaded ready-built - no source compilation), and fall back to CRAN source only if
    # that is unavailable or the binary will not load in this (conda) R. Either way build
    # with every core (Ncpus / MAKEFLAGS -j). The `requireNamespace` check inside the script
    # loads the namespace, so an ABI-mismatched binary is detected and the source build runs.
    # `|| true` + the verify below keep a harmless non-zero exit from killing setup under set -e.
    LF_R="$(mktemp --suffix=.R 2>/dev/null || echo "/tmp/lnckit_lncfinder.$$.R")"
    cat > "$LF_R" <<'RS'
os  <- tryCatch(readLines("/etc/os-release"), error=function(e) character(0))
val <- function(k){ v<-grep(paste0("^",k,"="),os,value=TRUE); if(length(v)) gsub('"','',sub(paste0("^",k,"="),"",v[1])) else "" }
id<-val("ID"); code<-val("VERSION_CODENAME"); vid<-val("VERSION_ID")
p3m <- ""
if (id %in% c("ubuntu","debian") && nzchar(code)) p3m <- sprintf("https://packagemanager.posit.co/cran/__linux__/%s/latest", code)
if (id %in% c("rhel","centos","rocky","almalinux","fedora")) { v<-substr(vid,1,1); if (v %in% c("7","8","9")) p3m <- sprintf("https://packagemanager.posit.co/cran/__linux__/rhel%s/latest", v) }
ncpus <- max(1L, tryCatch(parallel::detectCores(), error=function(e) 1L))
ua <- sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version$platform, R.version$arch, R.version$os))
do_install <- function(repo, binary){
  options(repos=c(CRAN=repo), HTTPUserAgent=if (binary) ua else sprintf("R (%s)", getRversion()))
  try(install.packages("LncFinder", dependencies=FALSE, Ncpus=ncpus), silent=TRUE)
}
if (nzchar(p3m)) do_install(p3m, TRUE)                                              # prebuilt binary: no compile
if (!requireNamespace("LncFinder", quietly=TRUE)) do_install("https://cloud.r-project.org", FALSE)  # source fallback
RS
    MAKEFLAGS="-j$(nproc 2>/dev/null || echo 2)" "$LF_RSCRIPT" --vanilla "$LF_R" || true
    rm -f "$LF_R"
    if "$LF_RSCRIPT" -e 'quit(status=!requireNamespace("LncFinder",quietly=TRUE))' 2>/dev/null; then
        echo "[ok]    LncFinder installed"
    else
        echo "[warn]  LncFinder did not install - the ensemble will just skip it"
    fi
fi

# ---- 4. FEELnc libdb compat fix -------------------------------------------
FEELNC_LIB="$ENV_ROOT/lnckit_feelnc/lib"
# `-print -quit` stops at the first match - avoids `... | head -1`, which under
# `pipefail` makes find exit non-zero (SIGPIPE) and trips `set -e`.
DBFILE_SO=$(find "$FEELNC_LIB/perl5" -name DB_File.so -print -quit 2>/dev/null || true)
if [ -n "$DBFILE_SO" ] && ldd "$DBFILE_SO" 2>/dev/null | grep -q "libdb-6.1.so => not found" \
   && [ ! -e "$FEELNC_LIB/libdb-6.1.so" ]; then
    echo "[fix]   FEELnc DB_File needs libdb-6.1.so - fetching a matching copy"
    TMP=$(mktemp -d)
    # Optional component: a failed fetch must not abort setup (FEELnc just won't run).
    if create_prefix -p "$TMP/libdb61" -c conda-forge libdb=6.1.26 >/dev/null 2>&1 \
       && [ -e "$TMP/libdb61/lib/libdb-6.1.so" ]; then
        cp "$TMP/libdb61/lib/libdb-6.1.so" "$FEELNC_LIB/libdb-6.1.so"
        echo "[fix]   installed $FEELNC_LIB/libdb-6.1.so"
    else
        echo "[warn]  could not fetch libdb-6.1.so - FEELnc DB_File may fail; use 'classifier: gffcompare'"
    fi
    rm -rf "$TMP"
else
    echo "[skip]  FEELnc libdb compat already satisfied"
fi

# ---- 4b. pre-install the pipeline TOOLS so ./lnckit.sh never stops to install ---
# The tools NOT baked into the conda env (SRA Toolkit, CNCI, PLEK, BLAST,
# Trimmomatic, extra aligners/assemblers) are otherwise installed on the FIRST
# lnckit.sh. We install them NOW instead, using the very same doctor()/installer the
# run uses. Three scopes:
#   default             -> what the DEFAULT pipeline needs (STANDARD mode) + SRA
#   LNCKIT_CONFIG=<yaml> -> exactly what YOUR pipeline.yaml selects (custom tools)
#   LNCKIT_PREINSTALL_ALL=1 -> every tool LncKit knows (all aligners/assemblers/...)
# Add specific extras with LNCKIT_PREINSTALL="blast trimmomatic".
# Skip the whole step with LNCKIT_SKIP_TOOLS=1 (tools then install on first run).
MANAGED_BIN="${LNCKIT_HOME:-${XDG_DATA_HOME:-$HOME/.local/share}}/lnckit/bin"
if [ "${LNCKIT_SKIP_TOOLS:-0}" = "1" ]; then
    echo "[skip]  tool pre-install (LNCKIT_SKIP_TOOLS=1) - tools install on first run"
else
    echo "[tools] pre-installing pipeline tools so the first ./lnckit.sh doesn't have to"
    # Put ALL env bins on PATH so doctor() sees tools already installed as present - this
    # includes the ISOLATED envs (lnckit_feelnc, lnckit_lncfinder), otherwise FEELnc and
    # LncFinder look 'missing' and get force-installed into the main lnckit env, where
    # FEELnc's old R/openssl/glibc deps make the solve fail. Only genuinely-missing tools
    # that belong in the lnckit env are installed here.
    PATH="$ENV_ROOT/lnckit/bin:$ENV_ROOT/lnckit_feelnc/bin:$ENV_ROOT/lnckit_lncfinder/bin:$MANAGED_BIN:$PATH" \
        PYTHONPATH="$HERE" \
        "$ENV_ROOT/lnckit/bin/python" - <<'PY' || \
        echo "[warn]  tool pre-install hit an error - the first ./lnckit.sh will retry it"
import os
from lnckit.core.config import load_config
from lnckit.core.installer import doctor, execute_plan, plan_install

if os.environ.get("LNCKIT_PREINSTALL_ALL") == "1":
    statuses = doctor(None)                       # every known tool
else:
    statuses = doctor(load_config(os.environ.get("LNCKIT_CONFIG") or None))

need = {s.tool for s in statuses if not s.present}
if os.environ.get("LNCKIT_SKIP_SRA") != "1":      # SRA is the default input path
    need |= {"prefetch", "fastq-dump"}
need |= set((os.environ.get("LNCKIT_PREINSTALL") or "").split())
# These live in their OWN conda envs (built by setup.sh's dedicated steps), NOT the main
# lnckit env - never try to conda-install them here (FEELnc's deps break the lnckit solve).
need -= {"feelnc", "lncfinder"}
need = sorted(need)

if not need:
    print("[tools] all pipeline tools already present")
else:
    print("[tools] installing:", ", ".join(need))
    # env="lnckit": conda packages (blast, trimmomatic, ...) land IN the lnckit env,
    # on lnckit.sh's PATH. Source/tarball tools (SRA, CNCI, PLEK) go to the managed bin.
    for line in execute_plan(plan_install(need), env="lnckit", yes=True):
        print("     ", line)
PY
fi

# ---- 4b-mqc. MultiQC optional deps (boto3 + pyarrow) ------------------------
# MultiQC 1.35 declares boto3 (S3 output) and pyarrow (Parquet) as requirements, but the
# conda build does not pull them. Two reasons the env can lack them: `conda env create`
# SKIPS an env that already exists (so updating envs/lnckit.yml alone never reaches an env
# built earlier), and a user may have created the env before they were added. Install them
# here explicitly and idempotently so the env is complete on BOTH fresh and existing
# installs - and so the CPAT `pip install` below stops reporting the unsatisfied requirement.
if [ "${LNCKIT_SKIP_TOOLS:-0}" != "1" ]; then
    if "$ENV_ROOT/lnckit/bin/python" -c 'import boto3, pyarrow' 2>/dev/null; then
        echo "[mqc]   boto3 + pyarrow already present"
    else
        echo "[mqc]   installing MultiQC optional deps (boto3 + pyarrow)"
        "$ENV_ROOT/lnckit/bin/python" -m pip install --quiet boto3 pyarrow \
            && echo "[mqc]   boto3 + pyarrow ready" \
            || echo "[warn]  could not install boto3/pyarrow (harmless - MultiQC still works)"
    fi
fi

# ---- 4b-cpat. Guarantee a Python-3 CPAT (v3) in the lnckit env --------------
# bioconda's solver sometimes lands the OLD CPAT 1.2.x, whose source is Python 2.
# Installed into the Python-3 lnckit env, its shebang runs it under python3 and it
# dies on 'print >>sys.stderr' (py2 syntax). We detect that broken state and force the
# maintained PyPI 'CPAT' (v3, Python 3 - same -x/-d/-g/-o CLI). Non-fatal: the coding
# ensemble just skips CPAT if this can't complete.
if [ "${LNCKIT_SKIP_TOOLS:-0}" != "1" ]; then
    CPAT_PY="$ENV_ROOT/lnckit/bin/cpat.py"
    broken=0
    if [ ! -f "$CPAT_PY" ]; then
        broken=1                                   # not installed at all
    elif grep -q 'print >>' "$CPAT_PY" 2>/dev/null; then
        broken=1                                   # Python-2 source -> crashes under py3
    fi
    if [ "$broken" = "1" ]; then
        echo "[cpat]  installing/repairing Python-3 CPAT (v3) via pip"
        "$ENV_ROOT/lnckit/bin/python" -m pip install --upgrade --force-reinstall CPAT \
            && echo "[cpat]  CPAT v3 ready" \
            || echo "[warn]  CPAT pip install failed - CPAT will be skipped in the ensemble"
    else
        echo "[cpat]  Python-3 CPAT already present"
    fi
fi

# ---- 4b-sra. Guarantee a MODERN SRA-Toolkit (>=3) is available ---------------
# Old sra-tools 2.x reject `prefetch --output-directory` and fail NCBI's TLS check, so
# accession download breaks. SRA Toolkit is a MANAGED-BIN static tarball (NCBI's official
# >=3 build), NOT a conda package: the bioconda sra-tools pins an old htslib that downgrades
# samtools in the shared env (see envs/lnckit.yml). So check where it REALLY is (managed bin,
# then PATH/env), and if it is somehow missing/old, re-fetch the static tarball with the SAME
# installer section 4b used - never `conda install` sra-tools into the env (that breaks
# samtools). Non-fatal: local-FASTQ runs do not need SRA at all.
if [ "${LNCKIT_SKIP_TOOLS:-0}" != "1" ]; then
    PF=""
    for c in "$MANAGED_BIN/prefetch" "$ENV_ROOT/lnckit/bin/prefetch"; do
        if [ -x "$c" ]; then PF="$c"; break; fi
    done
    if [ -z "$PF" ]; then PF="$(command -v prefetch 2>/dev/null || true)"; fi
    ver=""
    if [ -n "$PF" ]; then
        ver="$("$PF" --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)"
    fi
    major="${ver%%.*}"
    if [ -n "$major" ] && [ "$major" -ge 3 ] 2>/dev/null; then
        echo "[sra ]  SRA Toolkit $ver (>=3) OK"
    else
        echo "[sra ]  SRA Toolkit is ${ver:-missing} - (re)fetching NCBI static tarball (>=3)"
        PATH="$ENV_ROOT/lnckit/bin:$MANAGED_BIN:$PATH" PYTHONPATH="$HERE" \
            "$ENV_ROOT/lnckit/bin/python" - <<'PY' \
            && echo "[sra ]  SRA Toolkit ready" \
            || echo "[warn]  could not install SRA Toolkit - download reads manually or use local FASTQ"
from lnckit.core.installer import execute_plan, plan_install
for line in execute_plan(plan_install(["prefetch", "fastq-dump"]), env="lnckit", yes=True):
    print("     ", line)
PY
    fi
fi

# ---- 4c. PLEK smoke test: prove the Python-2 build actually RUNS -------------
# PLEK is a Python-2 tool whose SVM helpers compile at install; a broken py2 env or
# a failed compile can leave a PLEK that installs but errors the moment it runs. We
# run it once here on a tiny generated FASTA so that failure surfaces NOW, at setup,
# not mid-pipeline. Non-fatal: the coding ensemble simply skips PLEK if it can't run.
if [ "${LNCKIT_SKIP_TOOLS:-0}" != "1" ]; then
    PATH="$ENV_ROOT/lnckit/bin:$MANAGED_BIN:$PATH" PYTHONPATH="$HERE" \
        "$ENV_ROOT/lnckit/bin/python" - <<'PY' || true
import os, re, shutil, subprocess, tempfile
# First, heal any stale Python-2 launcher (e.g. a PLEK/CNCI shim written by an older
# install that points at an env path which no longer exists - the classic
# '.../.conda/envs/lnckit_py2/bin/python: No such file' error). This only rewrites the
# shim in place (no rebuild), so the smoke test below exercises the FIXED launcher.
try:
    from lnckit.core.installer import repair_py2_shims
    _fixed = repair_py2_shims()
    if _fixed:
        print("[fix ]  repaired stale Python-2 launcher(s): %s" % ", ".join(_fixed))
except Exception as _e:
    print("[note]  could not run shim repair (%s) - continuing" % _e)
plek = shutil.which("PLEK") or shutil.which("PLEK.py")
if not plek:
    print("[skip]  PLEK smoke test - PLEK not installed")
else:
    # The managed-bin PLEK is a shim: `exec python2 /.../PLEK.py "$@"`. Run from the
    # real PLEK dir so PLEK finds its bundled model/range files.
    plek_dir = None
    try:
        # the shim quotes the path (exec python2 "/.../PLEK.py" "$@"); strip quotes so the
        # captured dir isn't like '"/.../PLEK.1.2'.
        m = re.search(r'([\'"]?)(\S+PLEK\.py)\1', open(plek).read())
        if m:
            plek_dir = os.path.dirname(m.group(2).strip('"\''))
    except Exception:
        pass
    if plek_dir and not os.path.isdir(plek_dir):
        print("[skip]  PLEK smoke test - PLEK not built yet (builds on first ./lnckit.sh)")
        raise SystemExit(0)
    d = tempfile.mkdtemp()
    try:
        import random
        random.seed(1)
        fa = os.path.join(d, "test.fa")
        with open(fa, "w") as fh:
            for i in range(3):
                fh.write(">t%d\n%s\n" % (i, "".join(random.choice("ACGT")
                                                    for _ in range(600))))
        out = os.path.join(d, "plek_out")
        r = subprocess.run([plek, "-fasta", fa, "-out", out, "-thread", "1"],
                           cwd=plek_dir, capture_output=True, text=True, timeout=600)
        if os.path.exists(out) and os.path.getsize(out) > 0:
            print("[ok]    PLEK smoke test passed - ran on a tiny FASTA, produced output")
        else:
            print("[warn]  PLEK installed but did NOT produce output - the ensemble will "
                  "skip PLEK. Last error:")
            print("        " + ((r.stderr or r.stdout or "").strip().splitlines() or
                                 ["(none)"])[-1][:200])
    except Exception as e:
        print("[warn]  PLEK smoke test could not run (%s) - the ensemble will skip PLEK"
              % e)
    finally:
        shutil.rmtree(d, ignore_errors=True)
PY
fi

# ---- 4c-bis. CPC2 (vendored, offline) - compile libsvm + shim under lnckit_py2 --
# The conda cpc2 ships without working svm-scale/svm-predict, so LncKit uses the
# VENDORED CPC2_standalone in third_party/ (Python 2; needs biopython+numpy = the
# lnckit_py2 env). We compile its bundled libsvm and drop a CPC2.py shim in the managed
# bin that runs it under the py2 env. Skipped cleanly if the py2 env isn't built yet.
CPC2_DIR="$HERE/third_party/CPC2_standalone"
CPC2_SVM="$CPC2_DIR/libs/libsvm/libsvm-3.18/svm-scale"
if [ -f "$CPC2_DIR/bin/CPC2.py" ]; then
    if [ ! -x "$CPC2_SVM" ]; then
        echo "[cpc2]  compiling bundled libsvm"
        ( cd "$CPC2_DIR/libs/libsvm" && tar xzf libsvm-3.18.tar.gz \
          && make -C libsvm-3.18 >/dev/null 2>&1 ) \
          && echo "[cpc2]  libsvm built" \
          || echo "[warn]  CPC2 libsvm build failed (need make + g++); CPC2 will be skipped"
    fi
    PY2="$ENV_ROOT/lnckit_py2/bin/python"
    if [ -x "$CPC2_SVM" ] && [ -x "$PY2" ]; then
        mkdir -p "$MANAGED_BIN"
        printf '#!/usr/bin/env bash\nexec "%s" "%s" "$@"\n' \
            "$PY2" "$CPC2_DIR/bin/CPC2.py" > "$MANAGED_BIN/CPC2.py"
        chmod +x "$MANAGED_BIN/CPC2.py"
        echo "[cpc2]  installed CPC2.py shim (vendored, runs under lnckit_py2)"
        # The broken conda cpc2 (if an older env installed it) sits FIRST on lnckit.sh's PATH
        # and would shadow this shim. Remove its scripts from the lnckit env so the working
        # vendored CPC2 wins. (rm is instant; a full `conda remove` on HPC is slow.)
        for f in CPC2.py CPC2 CPC2_output_peptide.py; do
            rm -f "$ENV_ROOT/lnckit/bin/$f"
        done
        echo "[cpc2]  removed any shadowing conda CPC2 from the lnckit env"
    else
        echo "[note]  CPC2 shim deferred - needs the lnckit_py2 env (built with CNCI/PLEK)"
    fi
fi

# ---- 4c-ter. RNAplonc (vendored) - ensure its precompiled binary is executable -----
# The txCdsPredict binary ships in third_party/; git/tar can drop the +x bit. RNAplonc
# also needs java (in the lnckit env) and BioPerl (reused from lnckit_feelnc).
RNAPLONC_BIN="$HERE/third_party/RNAplonc/txCdsPredict"
if [ -f "$RNAPLONC_BIN" ]; then
    chmod +x "$RNAPLONC_BIN" 2>/dev/null && \
        echo "[rnaplonc] txCdsPredict marked executable (Weka/REPTree plant backend ready)"
fi

# ---- 4d. (opt-in) download the shared reference databases at install time ----
# UniProt (~280 MB) + Rfam.cm (~1.5 GB) are too big to bundle in the repo/tarball
# (GitHub caps files at 100 MB) and they update over time - so they are NOT shipped.
# But you CAN fetch them here, once, as part of install, with LNCKIT_WITH_DB=1. On an
# HPC, run setup.sh on the LOGIN node with this set; the compute job then needs no net.
if [ "${LNCKIT_WITH_DB:-0}" = "1" ]; then
    echo "[db]    downloading reference databases (UniProt + Rfam + Pfam) - large, one-time"
    # cmpress/hmmpress live in the lnckit env; expose them so Rfam.cm/Pfam-A.hmm press too.
    PATH="$ENV_ROOT/lnckit/bin:$MANAGED_BIN:$PATH" bash "$HERE/download_databases.sh" \
        || echo "[warn]  database download failed - run download_databases.sh later"
else
    echo "[db]    skipping database download (set LNCKIT_WITH_DB=1 to fetch UniProt +"
    echo "        Rfam + Pfam now, or run: bash download_databases.sh)"
fi

# ---- 5. record where the envs live, so lnckit.sh finds them (any manager) -----
{
    echo "LNCKIT_ENV_ROOT=$ENV_ROOT"
    echo "LNCKIT_PKG=$PKG"
} > "$HERE/.lnckit_env"

# ---- 6. VERIFY the tools actually installed (confirm the downloads worked) --
echo
echo "----------------------------------------------------------------"
echo "Verifying installed tools:"
MISSING=0
have() {  # $1=env  $2=binary (glob ok)  $3=label
    if ls "$ENV_ROOT/$1/bin/"$2 >/dev/null 2>&1; then
        echo "   [OK]      $3"
    else
        echo "   [MISSING] $3   (expected in $1 env)"
        MISSING=$((MISSING + 1))
    fi
}
# Some tools live in the MANAGED bin (vendored/source), not a conda env. have_managed
# checks the lnckit env first, then the managed bin - so a vendored tool (CPC2) isn't
# wrongly reported missing just because it isn't inside the conda env.
have_managed() {  # $1=env  $2=binary (glob ok)  $3=label
    if ls "$ENV_ROOT/$1/bin/"$2 >/dev/null 2>&1 || ls "$MANAGED_BIN/"$2 >/dev/null 2>&1; then
        echo "   [OK]      $3"
    else
        echo "   [MISSING] $3   (expected in $1 env or managed bin)"
        MISSING=$((MISSING + 1))
    fi
}
have lnckit 'hisat2'            "HISAT2      (alignment)"
have lnckit 'stringtie'         "StringTie   (assembly)"
have lnckit 'cpat*'             "CPAT        (coding potential)"
have_managed lnckit 'CPC2*'     "CPC2        (coding potential)"
have lnckit 'diamond'           "DIAMOND     (protein homology)"
have lnckit 'cmscan'            "Infernal    (Rfam removal)"
have lnckit 'getorf'            "getorf/EMBOSS (ORF)"
have lnckit 'samtools'          "samtools    (sort/index)"
have lnckit 'fastp'             "fastp       (trim/QC)"
have lnckit 'gffcompare'        "gffcompare  (compare)"
have lnckit 'gffread'           "gffread     (extract FASTA)"
if [ -d "$ENV_ROOT/lnckit_feelnc" ]; then
    have lnckit_feelnc 'FEELnc_classifier.pl' "FEELnc     (classification)"
else
    echo "   [skip]    FEELnc not installed - use 'classifier: gffcompare'"
fi
if [ -d "$ENV_ROOT/lnckit_lncfinder" ]; then
    have lnckit_lncfinder 'Rscript'  "LncFinder  (R, coding potential)"
fi

echo
if [ "$MISSING" -eq 0 ]; then
    echo "All tools installed OK. Environments under: $ENV_ROOT"
    echo "Note: CNCI & PLEK build on the first ./lnckit.sh (they are not conda packages)."
    echo "Now run:   ./lnckit.sh /path/to/your_project"
else
    echo "WARNING: $MISSING tool(s) missing above - the pipeline may skip those steps."
    echo "Re-run 'bash setup.sh' (downloads can fail on a flaky network), or check the"
    echo "output above for the failing package."
    exit 1
fi
