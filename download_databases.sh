#!/usr/bin/env bash
# download_databases.sh - fetch the shared reference databases ONCE, so every run (and
# every user pointing at the same LNCKIT_DB_DIR) reuses them without re-downloading.
#
# These files are large and update over time, so they are NOT bundled in the LncKit
# package/tarball - you download the current versions here instead. Run this on a machine
# (or HPC LOGIN node) that HAS internet; on an HPC whose compute nodes have no network,
# run it on the login node into shared storage and the compute job will just read it.
#
#   UniProt Swiss-Prot (protein DB for the homology filter)  ~90 MB gz  -> ~280 MB
#   Rfam.cm            (Infernal covariance models, ncRNA removal) ~350 MB gz -> ~1.5 GB
#   Pfam-A.fasta       (protein DB for the homology filter)  ~6.3 GB gz -> ~20 GB
#   COG 2020 proteins  (protein DB for the homology filter)  ~645 MB gz -> ~1.6 GB
#
# Usage (NO flags = download EVERYTHING; flags are ADDITIVE, to pick a subset):
#   bash download_databases.sh                       # DEFAULT: ALL (UniProt+Rfam+Pfam+COG) -> ./databases
#   LNCKIT_DB_DIR=/scratch/$USER/lnckit_db bash download_databases.sh   # shared location
#   bash download_databases.sh --uniprot --pfam --cog   # pick exactly what you want (subset)
#   bash download_databases.sh --homology            # all 3 homology DBs: UniProt + Pfam + COG
#   bash download_databases.sh --all                 # everything (same as no flags)
#   bash download_databases.sh --uniprot             # UniProt only          ~0.3 GB
#   bash download_databases.sh --rfam                # Rfam only             ~1.5 GB
#   bash download_databases.sh --cog                 # COG protein FASTA     ~1.6 GB
#   bash download_databases.sh --pfam                # Pfam-A protein FASTA  ~20 GB (BLASTX)
#   bash download_databases.sh --cpat=human          # a CPAT model (human|mouse|fly|zebrafish)
#   bash download_databases.sh --cpat                # ALL animal CPAT models
#     (no PLANT model in the CPAT repo - supply Plant_Hexamer.tsv + Plant.logit.RData)
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DB_DIR="${LNCKIT_DB_DIR:-$HERE/databases}"
# Parallel connections per file for aria2c (more = faster on per-connection-throttled mirrors
# like EBI). Override with LNCKIT_DL_CONN; LNCKIT_DL_TRIES sets how many from-scratch retries.
DL_CONN="${LNCKIT_DL_CONN:-16}"
# Flags are ADDITIVE: each one turns ON that database, so they combine freely, e.g.
#   bash download_databases.sh --uniprot --pfam --cog     # pick a SUBSET
# With NO flags, EVERYTHING is fetched: UniProt + Rfam + Pfam + COG (what the default config
# uses). This is large - Pfam-A alone is ~20 GB - so pass flags above to fetch only a subset
# if you do not need all of them. `--homology` grabs the three homology DBs (UniProt+Pfam+COG).
DO_UNIPROT=0; DO_RFAM=0; DO_COG=0; DO_PFAM_FASTA=0; CPAT_ORG=""; SELECTED=0
for a in "$@"; do
    SELECTED=1
    case "$a" in
        --uniprot)           DO_UNIPROT=1 ;;
        --rfam)              DO_RFAM=1 ;;
        --cog)               DO_COG=1 ;;
        --pfam|--pfam-fasta) DO_PFAM_FASTA=1 ;;        # Pfam-A protein FASTA (BLASTX)
        --homology)          DO_UNIPROT=1; DO_PFAM_FASTA=1; DO_COG=1 ;;  # all 3 homology DBs
        --all)               DO_UNIPROT=1; DO_RFAM=1; DO_PFAM_FASTA=1; DO_COG=1 ;;
        --cpat=*)  CPAT_ORG="${a#*=}" ;;
        --cpat)    CPAT_ORG="all" ;;
        -h|--help) sed -n '2,27p' "$0"; exit 0 ;;
        *) echo "unknown option: $a" >&2; exit 2 ;;
    esac
done
# No flags -> download EVERYTHING: UniProt + Rfam + Pfam + COG (the full default set).
if [ "$SELECTED" = 0 ]; then DO_UNIPROT=1; DO_RFAM=1; DO_PFAM_FASTA=1; DO_COG=1; fi
mkdir -p "$DB_DIR"
echo "Databases -> $DB_DIR"
if command -v aria2c >/dev/null 2>&1; then
    echo "[info] using aria2c ($DL_CONN connections/file) for faster downloads"
else
    echo "[info] aria2c not found - falling back to curl (single connection, slower)."
    echo "       For much faster (multi-connection) downloads, get aria2 onto PATH, e.g.:"
    echo "         micromamba install -n lnckit -c conda-forge aria2"
    echo "         micromamba run -n lnckit ./download_databases.sh   # run WITH the env so aria2c is found"
fi

# fetch <url> <final-path>: integrity-checked download that ALWAYS starts from the beginning.
#
# Policy (deliberately NOT resume-based): resuming a multi-GB partial over HTTP Range was
# unreliable on some mirrors - EBI/Pfam would stall for hours part-way through - so we never
# resume. Instead:
#   * REUSE a file only when it genuinely FINISHED: a <final>.done marker records its exact
#     current size, OR (for a .gz source) a complete <final>.part.gz is already on disk and
#     `gzip -t` proves it whole. Fully downloaded -> skipped.
#   * Anything else - a half-written partial from a broken/interrupted download, or a file
#     with no verifiable marker - is DELETED and downloaded again FROM BYTE 0. Half downloaded
#     -> overwritten and restarted from the start.
#   * For a .gz source the whole archive is verified with `gzip -t` BEFORE decompressing (this
#     catches a truncated download), then promoted to the final name with an atomic mv and its
#     verified byte-size recorded in <final>.done.
# So a connection that breaks mid-transfer means the next run starts that file over, not
# resumes it. The old final file is left untouched until a freshly verified copy replaces it.
#
# To bless a file you are SURE is complete (and skip re-downloading it), record its size
# yourself:   stat -c%s <final-path> > <final-path>.done

# dl <url> <out-file>: download FROM THE START (never resumes). Prefers aria2c, which opens
# several connections and splits the file across them - on per-connection-throttled mirrors
# like EBI this is several times faster than a single curl stream (DL_CONN connections; tune
# with LNCKIT_DL_CONN). Falls back to curl when aria2c is absent. Each attempt first DELETES
# any partial so the transfer always begins at byte 0; a few from-scratch retries ride out a
# flaky link, and on final give-up the partial is removed so a re-run simply starts clean.
dl() {
    local url="$1" out="$2" tries=0 max="${LNCKIT_DL_TRIES:-3}"
    while :; do
        tries=$((tries + 1))
        rm -f "$out" "$out.aria2"                     # ALWAYS start from byte 0 (no resume)
        if command -v aria2c >/dev/null 2>&1; then
            # -x/-s: open DL_CONN connections and split the file across them. --continue=false
            # + --allow-overwrite force a clean from-scratch download; --max-tries=1 because the
            # loop below owns retries (each one restarts cleanly).
            aria2c -x "$DL_CONN" -s "$DL_CONN" -k 1M \
                   --continue=false --allow-overwrite=true --auto-file-renaming=false \
                   --max-tries=1 --connect-timeout=30 \
                   --console-log-level=warn --summary-interval=0 \
                   -d "$(dirname "$out")" -o "$(basename "$out")" "$url" && return 0
        else
            # curl fallback: single connection, from byte 0 (NO -C -). --retry rides out brief
            # transient errors by re-issuing the request, which also restarts from the start.
            curl -fL --retry 2 --retry-delay 5 --connect-timeout 30 -o "$out" "$url" && return 0
        fi
        if [ "$tries" -ge "$max" ]; then
            rm -f "$out" "$out.aria2"                  # give up clean: leave no half file behind
            echo "[err ] $(basename "$out"): download failed after $max from-scratch attempt(s)." >&2
            echo "       Check internet/DNS reachability, then re-run - it restarts this file from the start." >&2
            return 1
        fi
        echo "[retry] $(basename "$out"): attempt $tries failed - deleting partial, restarting from byte 0 in 10s..." >&2
        sleep 10
    done
}

fetch() {
    local url="$1" out="$2" part="$2.part" done="$2.done"
    if [ -s "$out" ] && [ -f "$done" ] \
       && [ "$(stat -c%s "$out" 2>/dev/null)" = "$(cat "$done" 2>/dev/null)" ]; then
        echo "[skip] $(basename "$out") already complete (verified)"; return
    fi
    [ -e "$out" ] && echo "[warn] $(basename "$out") unverified (no matching .done) - re-downloading in full"
    rm -f "$done"                                    # drop any stale/unmatched marker
    echo "[get ] $(basename "$out")  <- $url"
    if [ "${url##*.}" = "gz" ]; then
        # gz source. Skip the download ONLY if a COMPLETE archive is already on disk (a prior
        # run finished downloading but was killed before decompress/promote); gzip -t proves it
        # whole, so re-fetching a huge file like Pfam-A would just waste a long transfer.
        # Otherwise a partial is a broken download -> delete it and fetch again from byte 0
        # (dl does the delete), then verify. No resume is ever attempted.
        if [ -s "$part.gz" ] && gzip -t "$part.gz" 2>/dev/null; then
            echo "[have] $(basename "$out"): complete archive already downloaded - skipping download"
        else
            [ -e "$part.gz" ] && echo "[fresh] $(basename "$out"): half-downloaded partial found - deleting it and downloading from the start" >&2
            dl "$url" "$part.gz" || exit 1               # downloads from byte 0 (no resume)
            if ! gzip -t "$part.gz" 2>/dev/null; then
                rm -f "$part.gz"                         # leave nothing half-written
                echo "[err ] $(basename "$out"): incomplete/corrupt after download - re-run to restart from the start." >&2
                exit 1
            fi
        fi
        gunzip -c "$part.gz" > "$part"
        rm -f "$part.gz"
    else
        # Raw (non-gz) file, e.g. a CPAT model: download from the start; dl deletes any partial
        # first and retries from byte 0, so a broken transfer is never resumed.
        dl "$url" "$part" || { echo "[err ] $(basename "$out"): download failed - re-run to restart from the start." >&2; exit 1; }
    fi
    [ -s "$part" ] || { echo "[err ] $(basename "$out"): empty after download." >&2; exit 1; }
    mv -f "$part" "$out"
    stat -c%s "$out" > "$done"                           # record verified size for skip checks
    echo "[ok  ] $out ($(stat -c%s "$out") bytes, verified)"
}

if [ "$DO_UNIPROT" = 1 ]; then
    fetch "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz" \
          "$DB_DIR/uniprot_sprot.fasta"
fi

# COG 2020 protein sequences (NCBI) for the diamond blastx homology filter (~645 MB gz).
# Add to reference.protein_db so a transcript hitting a COG protein is dropped as coding.
if [ "$DO_COG" = 1 ]; then
    fetch "https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/cog-20.fa.gz" "$DB_DIR/cog-20.fa"
    echo "COG in $DB_DIR/cog-20.fa - add it to reference.protein_db (a list) in your config."
fi

# Pfam-A protein SEQUENCES (EBI) for diamond blastx (~6.3 GB gz -> ~20 GB). Add to
# reference.protein_db so a transcript hitting a Pfam-A sequence is dropped as coding.
# Large: only fetch if you want Pfam in the homology filter.
if [ "$DO_PFAM_FASTA" = 1 ]; then
    fetch "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.fasta.gz" \
          "$DB_DIR/Pfam-A.fasta"
    echo "Pfam-A.fasta in $DB_DIR - set reference.pfam_protein to it (diamond blastx)."
fi

# CPAT prebuilt models (Liguo Wang's repo). Case matches the repo EXACTLY: fly/zebrafish
# hexamers are lowercase, everything else capitalized - a wrong case 404s. There is NO
# plant model here; supply Plant_Hexamer.tsv + Plant.logit.RData yourself for plants.
if [ -n "$CPAT_ORG" ]; then
    CPAT_REPO="https://raw.githubusercontent.com/liguowang/cpat/master/prebuilt_models"
    cpat_one() {   # $1 = hexamer filename, $2 = logit filename
        fetch "$CPAT_REPO/$1" "$DB_DIR/$1"
        fetch "$CPAT_REPO/$2" "$DB_DIR/$2"
    }
    case "$CPAT_ORG" in
        human|all)      cpat_one Human_Hexamer.tsv     Human_logitModel.RData ;;&
        mouse|all)      cpat_one Mouse_Hexamer.tsv     Mouse_logitModel.RData ;;&
        fly|all)        cpat_one fly_Hexamer.tsv       Fly_logitModel.RData ;;&
        zebrafish|all)  cpat_one zebrafish_Hexamer.tsv Zebrafish_logitModel.RData ;;&
        human|mouse|fly|zebrafish|all) : ;;   # matched at least one above
        *) echo "unknown --cpat organism '$CPAT_ORG' (human|mouse|fly|zebrafish|all)"; exit 2 ;;
    esac
    echo "CPAT model(s) for '$CPAT_ORG' in $DB_DIR - set organism accordingly (or reference.cpat_*)."
fi

if [ "$DO_RFAM" = 1 ]; then
    fetch "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz" "$DB_DIR/Rfam.cm"
    # Infernal needs the .cm "pressed" before cmscan can use it (build the binary indexes).
    if command -v cmpress >/dev/null 2>&1 && [ ! -s "$DB_DIR/Rfam.cm.i1i" ]; then
        echo "[press] cmpress Rfam.cm"
        if cmpress -F "$DB_DIR/Rfam.cm" >/dev/null; then
            echo "[ok  ] Rfam.cm pressed"
        else
            echo "[warn] cmpress failed - run 'cmpress $DB_DIR/Rfam.cm' manually before the housekeeping step"
        fi
    else
        [ -s "$DB_DIR/Rfam.cm.i1i" ] || echo "[note] run 'cmpress $DB_DIR/Rfam.cm' once " \
            "(cmpress is in the lnckit env) before the housekeeping step"
    fi
fi
# NOTE: Pfam is used as a BLASTX homology database (--pfam-fasta -> Pfam-A.fasta), not as an
# hmmscan coding backend, so there is no Pfam-A.hmm download here.

echo "Done. Point runs at these with:  export LNCKIT_DB_DIR=$DB_DIR"
