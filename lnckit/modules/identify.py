"""The three evidence stages that turn non-coding candidates into lncRNAs.

Each function runs one external tool over the *current* candidate FASTA and returns
a FilterResult, or None when the reference data that stage needs is not configured
(so a run without, say, an Rfam database still completes - loudly, with a warning,
rather than silently pretending the filter was applied).

    coding_potential : could this transcript encode a protein?  (CPAT / CPC2 / PLEK)
    protein_homology : does it align to a known protein?        (DIAMOND / BLASTX)
    housekeeping     : is it a known structural ncRNA?          (Infernal+Rfam / BLASTN)

The third is the one Plant-LncPipe omits, and it is not redundant with the first
two: rRNA, tRNA, snoRNA and miRNA precursors are *genuinely non-coding*, so they
pass coding-potential and protein-homology filtering unscathed and end up in the
final set unless they are explicitly removed here.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional, Set

from lnckit.core.context import RunContext
from lnckit.core.logging_utils import get_logger
from lnckit.filters import coding_potential as cp
from lnckit.filters.base import FilterResult
from lnckit.filters.housekeeping import (
    HousekeepingConfig, filter_housekeeping, parse_blast_hits, parse_cmscan_tblout)
from lnckit.filters.protein_homology import filter_protein_hits
from lnckit.utils.fasta import subset_fasta
from lnckit.wrappers import get_wrapper

log = get_logger("lnckit.identify")


def candidate_fasta(ctx: RunContext, tx_fasta: Path, candidates: Set[str],
                    tag: str) -> Path:
    """Write the surviving candidates to their own FASTA (keeps each tool's input small)."""
    out = ctx.dirs["filtered_transcripts"] / f"candidates_{tag}.fasta"
    return subset_fasta(tx_fasta, candidates, out)


def _empty_result(name: str) -> FilterResult:
    """No candidates left: return an empty pass-through so the stage is a clean 0 -> 0
    instead of running a tool on an empty FASTA (cmscan/DIAMOND abort on empty input)."""
    log.info("%s: no candidates remain - nothing to test", name)
    return FilterResult(name, set(), set())


def _is_plant(ctx: RunContext) -> bool:
    """True when the configured organism is a plant (drives model defaults)."""
    return str(ctx.config.get("organism", "other")).strip().lower() in ("plant", "plants", "pl")


# Per-organism CPAT models. Explicit reference.cpat_hexamer/cpat_logit always win; when
# they are unset, the model is chosen from `organism`. Animal models are CPAT's prebuilt
# ones (auto-downloaded when the network allows); the PLANT model is the plant-retrained
# CPAT model (Plant_Hexamer.tsv + Plant.logit.RData) - place it in $LNCKIT_DB_DIR (no
# single canonical URL). Tuple: (hexamer_name, logit_name, hexamer_url, logit_url, cutoff).
_CPAT_REPO = "https://raw.githubusercontent.com/liguowang/cpat/master/prebuilt_models/"
_CPAT_MODELS = {
    "plant":     ("Plant_Hexamer.tsv", "Plant.logit.RData", None, None, 0.46),
    "human":     ("Human_Hexamer.tsv", "Human_logitModel.RData",
                  _CPAT_REPO + "Human_Hexamer.tsv", _CPAT_REPO + "Human_logitModel.RData", 0.364),
    "mouse":     ("Mouse_Hexamer.tsv", "Mouse_logitModel.RData",
                  _CPAT_REPO + "Mouse_Hexamer.tsv", _CPAT_REPO + "Mouse_logitModel.RData", 0.44),
    # NB the CPAT repo mixes case: fly/zebrafish HEXAMERs are lowercase, logit models
    # capitalized - filenames below match the repo exactly (a wrong case 404s).
    "fly":       ("fly_Hexamer.tsv", "Fly_logitModel.RData",
                  _CPAT_REPO + "fly_Hexamer.tsv", _CPAT_REPO + "Fly_logitModel.RData", 0.39),
    "zebrafish": ("zebrafish_Hexamer.tsv", "Zebrafish_logitModel.RData",
                  _CPAT_REPO + "zebrafish_Hexamer.tsv", _CPAT_REPO + "Zebrafish_logitModel.RData", 0.381),
}
# Map the free-text `organism` value to a model key.
_CPAT_ORG = {
    "plant": "plant", "plants": "plant", "pl": "plant",
    "human": "human", "hsapiens": "human", "vertebrate": "human", "animal": "human",
    "mouse": "mouse", "mmusculus": "mouse",
    "fly": "fly", "drosophila": "fly",
    "zebrafish": "zebrafish", "danio": "zebrafish",
}

# LncFinder ships built-in SVM models only for human and mouse (lnc_finder()'s own default is
# "human"). Map our `organism` aliases to the exact built-in name, used for BOTH svm.model and
# frequencies.file. Organisms with no LncFinder built-in (fly/zebrafish/other) are left to the
# wrapper's default; plant has its own retrained model handled below.
_LNCFINDER_BUILTIN = {
    "human": "human", "hsapiens": "human", "vertebrate": "human", "animal": "human",
    "mouse": "mouse", "mmusculus": "mouse",
}


def _db_dir(ctx: RunContext) -> Path:
    """Where shared reference models live (mirrors lnckit.sh's LNCKIT_DB_DIR)."""
    import os
    return Path(os.environ.get("LNCKIT_DB_DIR", "databases"))


def _download_if_missing(dest: Path, url: Optional[str]) -> None:
    """Best-effort fetch of a model file; stays silent-safe on a no-network node."""
    if dest.exists() or not url:
        return
    import urllib.request
    try:
        dest.parent.mkdir(parents=True, exist_ok=True)
        log.info("downloading CPAT model %s", url)
        urllib.request.urlretrieve(url, dest)          # noqa: S310 - curated CPAT URL
    except Exception as e:                              # noqa: BLE001 - offline nodes etc.
        log.warning("could not download CPAT model %s (%s) - place the file manually", url, e)


def _resolve_lncfinder_model(ctx: RunContext):
    """Return (model, frequencies) for LncFinder, chosen from `organism` - mirroring CPAT.

    Precedence: explicit coding_potential.lncfinder.model wins; else, for a PLANT organism,
    auto-use the bundled plant SVM model (Plant_model.rda) with frequencies recomputed from
    the bundled training FASTAs (make_frequencies) - exactly the plant LncFinder recipe. For
    any other organism, return (None, None) so the wrapper uses LncFinder's built-in model.
    So selecting `organism: plant` "just works" like CPAT's plant model, with no extra config.
    """
    lf = ctx.config.get("coding_potential.lncfinder", {}) or {}
    model = lf.get("model")
    if model:                                  # explicit config always wins
        return model, (lf.get("frequencies") or model)
    if not _is_plant(ctx):
        # Non-plant: use LncFinder's matching built-in model when it ships one for this
        # organism (human -> "human", mouse -> "mouse"); otherwise fall back to the wrapper's
        # default. The built-in name drives both svm.model and frequencies.file.
        org = str(ctx.config.get("organism", "other")).strip().lower()
        builtin = _LNCFINDER_BUILTIN.get(org)
        if builtin:
            log.info("LncFinder model for organism '%s': built-in '%s' model", org, builtin)
            return builtin, builtin
        return None, None                      # no built-in for this organism -> wrapper default
    # Plant: find Plant_model.rda + the training FASTAs in $LNCKIT_DB_DIR or bundled databases/
    bundled = Path(__file__).resolve().parents[2] / "databases"
    for d in (_db_dir(ctx), bundled):
        pm, trm, trl = d / "Plant_model.rda", d / "training_mRNA.fasta", d / "training_lncRNA.fasta"
        if pm.exists() and trm.exists() and trl.exists():
            log.info("LncFinder model for organism 'plant': %s + make_frequencies(%s, %s)",
                     pm.name, trm.name, trl.name)
            return str(pm), f"make:{trm},{trl}"
    log.warning("organism is 'plant' but the plant LncFinder model/training files "
                "(Plant_model.rda + training_mRNA.fasta + training_lncRNA.fasta) were not "
                "found in $LNCKIT_DB_DIR or databases/ - falling back to the built-in model.")
    return None, None


def _resolve_cpat_model(ctx: RunContext):
    """Return (hexamer, logit, cutoff) for CPAT, or (None, None, None) if unavailable.

    Precedence: explicit reference.cpat_hexamer/cpat_logit  ->  the model for `organism`
    (files under $LNCKIT_DB_DIR, animal ones auto-downloaded when reachable).
    """
    hexamer = ctx.config.get("reference.cpat_hexamer")
    logit = ctx.config.get("reference.cpat_logit")
    # A user-set coding_cutoff (a NUMBER) overrides the per-organism default; null/None means
    # "use the per-organism default", so read it as None-when-unset and apply the fallback below.
    user_cut = ctx.config.get("coding_potential.cpat.coding_cutoff")
    if hexamer and logit:
        return hexamer, logit, (user_cut if user_cut is not None else 0.46)
    org = _CPAT_ORG.get(str(ctx.config.get("organism", "other")).strip().lower())
    if not org:
        return None, None, None
    hx_name, lg_name, hx_url, lg_url, cutoff = _CPAT_MODELS[org]
    # Look for the model files in $LNCKIT_DB_DIR AND in the package's own bundled
    # databases/ dir. The PLANT model ships with LncKit (in databases/), so 'organism:
    # plant' works out of the box; animal models auto-download into $LNCKIT_DB_DIR.
    bundled = Path(__file__).resolve().parents[2] / "databases"
    search_dirs = [_db_dir(ctx), bundled]
    hx = next((d / hx_name for d in search_dirs if (d / hx_name).exists()), None)
    lg = next((d / lg_name for d in search_dirs if (d / lg_name).exists()), None)
    if not (hx and lg):                       # not found locally -> try downloading (animals)
        hx = _db_dir(ctx) / hx_name
        lg = _db_dir(ctx) / lg_name
        _download_if_missing(hx, hx_url)
        _download_if_missing(lg, lg_url)
    if hx and lg and hx.exists() and lg.exists():
        log.info("CPAT model for organism '%s': %s + %s", org, hx.name, lg.name)
        # a user-set cutoff (number) wins; null/None -> the per-organism default `cutoff`
        return str(hx), str(lg), (user_cut if user_cut is not None else cutoff)
    return None, None, None


# --------------------------------------------------------------------------- #
# 1. Coding potential
# --------------------------------------------------------------------------- #
def run_coding_potential(ctx: RunContext, tx_fasta: Path,
                         candidates: Set[str]) -> Optional[FilterResult]:
    """Run the configured coding-potential backend(s) and keep the non-coding set.

    `tools.coding_potential` may be a single tool or a **list**. A list reproduces
    the Plant-LncPipe ensemble: each tool is run on the same candidates and the
    non-coding sets are INTERSECTED, so a transcript survives only if *every* model
    calls it non-coding. Backends that are skipped (missing reference data) drop out
    of the intersection rather than vetoing everything.
    """
    if not candidates:
        return _empty_result("coding_potential")
    backends = ctx.config.tools.get("coding_potential", "cpat")
    # RNAplonc is a PLANT-trained classifier. It stays in the default ensemble, but it only
    # actually runs for a plant organism; for anything else it is set aside (not used), so the
    # "default = all tools" selection is safe for non-plant data. This applies no matter how
    # rnaplonc got into the list (default.yaml, pipeline.yaml or an explicit --coding).
    if not _is_plant(ctx):
        org = ctx.config.get("organism", "other")
        if isinstance(backends, (list, tuple)):
            if any(str(b).strip().lower() == "rnaplonc" for b in backends):
                backends = [b for b in backends if str(b).strip().lower() != "rnaplonc"]
                log.info("rnaplonc is plant-only - set aside for organism '%s' "
                         "(not used in the coding-potential ensemble)", org)
        elif str(backends).strip().lower() == "rnaplonc":
            log.warning("rnaplonc is plant-only and the organism is '%s', and no other "
                        "coding-potential tool is configured - skipping the coding step", org)
            return None
    if isinstance(backends, (list, tuple)):
        results = []
        for b in backends:
            r = _safe_coding_backend(ctx, b, tx_fasta, candidates)
            if r is not None:
                results.append((b, r))
        if not results:
            return None
        return _combine_coding(ctx, results, candidates)
    return _safe_coding_backend(ctx, backends, tx_fasta, candidates)


def _combine_coding(ctx: RunContext, results, candidates: Set[str]) -> FilterResult:
    """Combine several coding-potential calls into one non-coding set.

    `coding_potential.consensus` chooses HOW the votes combine (default "all"):
      * "all"       -> unanimous: keep only transcripts EVERY tool called non-coding
                       (strict intersection; the safest, fewest false lncRNAs).
      * "majority"  -> keep transcripts a strict majority (> half) called non-coding,
                       so one dissenting tool cannot veto - e.g. 5 of 6 non-coding keeps it.
      * <integer N> -> keep transcripts at least N tools called non-coding.

    Only tools that actually ran count toward the vote (a skipped tool is neither a
    yes nor a no), so the denominator is the number of successful backends.
    """
    from collections import Counter
    n = len(results)
    votes: Counter = Counter()
    for _b, r in results:
        for t in (candidates & r.kept):          # this tool called t non-coding
            votes[t] += 1

    mode = str(ctx.config.get("coding_potential.consensus", "all")).strip().lower()
    if mode in ("majority", "consensus", "vote"):
        threshold = n // 2 + 1                    # strict majority
        rule = f"non-coding by majority (>= {threshold} of {n})"
    elif mode.isdigit():
        threshold = max(1, min(int(mode), n))
        rule = f"non-coding by >= {threshold} of {n}"
    else:                                         # "all" (default): unanimous
        threshold = n
        rule = f"non-coding by all {n}"

    kept = {t for t in candidates if votes[t] >= threshold}
    names = "+".join(b for b, _ in results)
    log.info("coding-potential ensemble [%s], %s: %d/%d kept",
             names, rule, len(kept), len(candidates))
    return FilterResult("coding_potential", kept, candidates - kept)


def _safe_coding_backend(ctx: RunContext, backend: str, tx_fasta: Path,
                         candidates: Set[str]) -> Optional[FilterResult]:
    """Run one coding-potential backend, but never let its failure kill the pipeline.

    A tool that is missing, un-installed, or errors out (e.g. LncFinder's R package
    absent, PLEK/CNCI Python-2 hiccups) is logged and dropped - in an ensemble it
    falls out of the intersection; as the sole backend the stage is simply skipped.
    This keeps a whole run (and its saved outputs) from being lost to one flaky tool.
    """
    try:
        return _run_one_coding_backend(ctx, backend, tx_fasta, candidates)
    except Exception as e:
        log.warning("coding-potential tool '%s' failed (%s) - skipping it. Install/fix "
                    "it, or drop it from tools.coding_potential.", backend, e)
        return None


# --------------------------------------------------------------------------- #
# RNAplonc - plant-trained REPTree (Weka) coding-potential backend (vendored)
# --------------------------------------------------------------------------- #
def _ensure_cpc2_libsvm() -> None:
    """Compile the vendored CPC2's bundled libsvm if svm-scale/svm-predict are missing.

    CPC2 shells out to libsvm's svm-scale/svm-predict at
    third_party/CPC2_standalone/libs/libsvm/libsvm-3.18/. Those are C programs built from
    a bundled tarball; if setup.sh didn't build them (e.g. re-run needed, or no compiler at
    install time) CPC2 fails with 'No excutable svm-scale on CPC2 path!'. We build them
    once here, so CPC2 self-heals. Best-effort: needs make + a C++ compiler; if that's not
    available, CPC2 is skipped (the ensemble continues) with a clear message.
    """
    import subprocess
    libsvm = (Path(__file__).resolve().parents[2]
              / "third_party" / "CPC2_standalone" / "libs" / "libsvm")
    if (libsvm / "libsvm-3.18" / "svm-scale").exists():
        return
    tarball = libsvm / "libsvm-3.18.tar.gz"
    if not tarball.exists():
        return
    log.info("CPC2: building bundled libsvm (svm-scale/svm-predict) - one-time")
    try:
        subprocess.run(["tar", "xzf", str(tarball)], cwd=str(libsvm), check=True)
        subprocess.run(["make", "-C", str(libsvm / "libsvm-3.18")], check=True,
                       capture_output=True)
    except Exception as e:                             # noqa: BLE001 - no compiler etc.
        log.warning("CPC2: could not build libsvm (%s). CPC2 will be skipped - install "
                    "make + g++ and re-run 'bash setup.sh', or drop cpc2 from "
                    "tools.coding_potential.", e)


def _rnaplonc_dir() -> Path:
    """Path to the vendored RNAplonc tool (third_party/RNAplonc at the repo root)."""
    return Path(__file__).resolve().parents[2] / "third_party" / "RNAplonc"


def _bioperl_perl():
    """A perl that has BioPerl (feature_extraction.pl needs Bio::DB::Fasta), plus its
    PERL5LIB. Reuse the isolated lnckit_feelnc env (it ships BioPerl); fall back to a
    plain `perl` on PATH. Returns (perl_exe, perl5lib_or_None) or (None, None)."""
    import shutil
    import sys
    feelnc = Path(sys.prefix).parent / "lnckit_feelnc"
    perl = feelnc / "bin" / "perl"
    if perl.exists():
        mods = sorted(feelnc.glob("lib/perl5/site_perl/*/Bio/SeqIO.pm"))
        return str(perl), (str(mods[0].parents[1]) if mods else None)
    p = shutil.which("perl")
    return (p, None) if p else (None, None)


def _run_rnaplonc(ctx: RunContext, fasta: Path, outdir: Path,
                  candidates: Set[str]) -> Optional[FilterResult]:
    """Run the vendored RNAplonc chain (plant REPTree/Weka) and label transcripts.

    Chain: txCdsPredict -> feature_extraction.pl (BioPerl) -> weka REPTree -l model -p 0.
    Experimental: any missing piece (java, BioPerl) makes it skip gracefully, and the
    ensemble simply continues without it.
    """
    import shutil
    from lnckit.utils.shell import run as run_tool
    from lnckit.filters import coding_potential as cp

    tool_dir = _rnaplonc_dir()
    txcds = tool_dir / "txCdsPredict"
    weka = tool_dir / "weka.jar"
    model = tool_dir / "RNAplonc.model"
    feat = tool_dir / "feature_extraction.pl"
    if not all(p.exists() for p in (txcds, weka, model, feat)):
        log.warning("coding-potential 'rnaplonc' skipped: vendored files missing in %s",
                    tool_dir)
        return None
    if not shutil.which("java"):
        log.warning("coding-potential 'rnaplonc' skipped: java not found "
                    "(install openjdk, e.g. add it to the lnckit env).")
        return None
    perl, perl5lib = _bioperl_perl()
    if not perl:
        log.warning("coding-potential 'rnaplonc' skipped: no perl available.")
        return None

    cds = outdir / "rnaplonc.cds"
    arff = outdir / "rnaplonc.arff"
    weka_out = outdir / "rnaplonc.weka.txt"
    try:
        txcds.chmod(txcds.stat().st_mode | 0o111)         # ensure the binary is executable
        rlog = ctx.logs / "rnaplonc.log"
        run_tool([str(txcds), str(fasta), str(cds)], tool="txCdsPredict",
                 log_file=rlog)
        run_tool([perl, str(feat), str(fasta), str(cds)], tool="rnaplonc-features",
                 log_file=rlog, stdout_path=arff,
                 env={"PERL5LIB": perl5lib} if perl5lib else None)
        run_tool(["java", "-cp", str(weka),
                  "weka.classifiers.trees.REPTree", "-l", str(model),
                  "-T", str(arff), "-p", "0"], tool="rnaplonc-weka",
                 log_file=rlog, stdout_path=weka_out)
    except Exception as e:                                 # noqa: BLE001 - experimental tool
        log.warning("coding-potential 'rnaplonc' failed (%s) - skipping it.", e)
        return None
    return cp.filter_by_label(cp.parse_rnaplonc(cds, weka_out, candidates),
                              candidates, "rnaplonc")


def _run_one_coding_backend(ctx: RunContext, backend: str, tx_fasta: Path,
                            candidates: Set[str]) -> Optional[FilterResult]:
    outdir = ctx.dirs["coding_prediction"]
    # Ensure the output dir exists - the coding tools (CPAT/CPC2/PLEK/...) write their
    # intermediate files here, and it may have been deleted between runs (e.g. a manual
    # 'rm -rf coding_prediction' to force a re-run), which crashed the tools with
    # 'No such file or directory'. (CNCI still makes its OWN subdir, so this only creates
    # the parent, which is safe.)
    outdir.mkdir(parents=True, exist_ok=True)
    fasta = candidate_fasta(ctx, tx_fasta, candidates, f"coding_{backend}")
    threads, logs = ctx.config.threads, ctx.logs

    if backend == "cpat":
        # Model is chosen from `organism` (plant -> Plant model, human/mouse/fly/zebrafish
        # -> that animal model), unless reference.cpat_hexamer/cpat_logit are set explicitly.
        hexamer, logit, cutoff = _resolve_cpat_model(ctx)
        if not (hexamer and logit):
            org = str(ctx.config.get("organism", "other")).strip().lower()
            log.warning(
                "coding-potential skipped: no CPAT model for organism '%s'. Set "
                "reference.cpat_hexamer + reference.cpat_logit (e.g. the plant model "
                "Plant_Hexamer.tsv + Plant.logit.RData), or place the model in "
                "$LNCKIT_DB_DIR, or set organism to plant/human/mouse/fly/zebrafish.", org)
            return None
        prefix = outdir / "cpat"
        # Optional explicit CPAT command (e.g. a CPAT 1.2.4 shim in the py2 env, to exactly
        # reproduce Plant-LncPipe). Default: whatever cpat.py/cpat is on PATH (v3).
        cpat_exe = ctx.config.get("coding_potential.cpat.executable") or None
        # Optional extra ORF-length gate (the ricebean paper used ORF_size < 100 on top of
        # coding_prob < 0.46). null -> coding_prob only.
        max_orf = ctx.config.get("coding_potential.cpat.max_orf_size")
        get_wrapper("coding_potential", "cpat")(threads=threads, log_dir=logs).run(
            fasta=fasta, out_prefix=prefix, hexamer=hexamer, logit_model=logit,
            executable=cpat_exe)
        return cp.filter_cpat(prefix, candidates, coding_cutoff=cutoff,
                              max_orf_size=float(max_orf) if max_orf is not None else None)

    if backend == "cpc2":
        _ensure_cpc2_libsvm()          # build the vendored libsvm if svm-scale is missing
        prefix = outdir / "cpc2"
        get_wrapper("coding_potential", "cpc2")(threads=threads, log_dir=logs).run(
            fasta=fasta, out_prefix=prefix)
        return cp.filter_by_label(cp.parse_cpc2(Path(f"{prefix}.txt")), candidates, "cpc2")

    if backend == "cnci":
        # CNCI creates its output directory itself and errors if it already exists,
        # so we pass the path but must NOT pre-create it. Model: explicit setting wins,
        # else pick from organism - CNCI has only 'pl' (plant) and 've' (vertebrate).
        out_dir = outdir / "cnci"
        model = ctx.config.get("coding_potential.cnci.model") or (
            "pl" if _is_plant(ctx) else "ve")
        get_wrapper("coding_potential", "cnci")(threads=threads, log_dir=logs).run(
            fasta=fasta, out_dir=out_dir, model=model)
        return cp.filter_by_label(cp.parse_cnci(out_dir / "CNCI.index"), candidates, "cnci")

    if backend == "plek":
        out_file = outdir / "plek.txt"
        get_wrapper("coding_potential", "plek")(threads=threads, log_dir=logs).run(
            fasta=fasta, out_file=out_file)
        return cp.filter_by_label(cp.parse_plek(out_file), candidates, "plek")

    if backend == "lncfinder":
        # Model chosen from `organism` (plant -> bundled Plant_model.rda + make_frequencies),
        # exactly like CPAT; non-plant uses LncFinder's built-in model. Explicit config wins.
        model, freq = _resolve_lncfinder_model(ctx)
        if not model:
            log.warning("LncFinder is using a built-in model (human/mouse). For a non-plant "
                        "organism set coding_potential.lncfinder.model (+ .frequencies) to a "
                        "model trained on it for best accuracy.")
        out_tsv = outdir / "lncfinder.tsv"
        get_wrapper("coding_potential", "lncfinder")(threads=threads, log_dir=logs).run(
            fasta=fasta, out_tsv=out_tsv, model=(model or "mouse"), frequencies=freq)
        return cp.filter_by_label(cp.parse_lncfinder(out_tsv), candidates, "lncfinder")

    if backend == "rnaplonc":
        return _run_rnaplonc(ctx, fasta, outdir, candidates)

    # NOTE: Pfam is intentionally NOT a coding-potential (hmmscan) backend. Pfam domains are
    # detected as a diamond BLASTX homology filter instead - add Pfam-A.fasta to
    # reference.protein_db (or reference.pfam_protein). See run_protein_filter below.

    log.warning("no coding-potential backend registered for '%s' - stage skipped", backend)
    return None


# --------------------------------------------------------------------------- #
# 2. Protein homology
# --------------------------------------------------------------------------- #
def _protein_databases(ctx: RunContext):
    """The protein FASTAs to search against, as (label, path). You can set either or
    both - a transcript is dropped if it matches a protein in ANY of them:
      * reference.protein_db     -> a broad DB (UniProt/SwissProt); may be a list
      * reference.genome_protein -> YOUR species' own annotated proteins
    """
    dbs = []
    pd = ctx.config.get("reference.protein_db")
    if isinstance(pd, (list, tuple)):
        dbs += [(f"uniprot{i+1}", p) for i, p in enumerate(pd) if p]
    elif pd:
        dbs.append(("uniprot", pd))
    # Pfam-A protein sequences (a FASTA, NOT the .hmm) - searched AFTER UniProt with the
    # SAME evalue/min_pident, removing transcripts that hit a Pfam family sequence.
    pf = ctx.config.get("reference.pfam_protein")
    if pf:
        dbs.append(("pfam", pf))
    gp = ctx.config.get("reference.genome_protein")
    if gp:
        dbs.append(("genome", gp))
    # Auto-include COG / Pfam-A protein FASTAs if they're sitting in the DB dir (downloaded
    # via download_databases.sh --cog / --pfam-fasta) but not already configured - so simply
    # downloading them adds them to the homology filter with no config edit. UniProt+Pfam+COG
    # then all run (each with pident>=min_pident & e<=evalue). Explicit config is still honored.
    configured = {str(Path(p)) for _, p in dbs if p}
    dbdir = _db_dir(ctx)
    for lbl, fname in (("cog", "cog-20.fa"), ("pfam", "Pfam-A.fasta")):
        f = dbdir / fname
        if f.exists() and str(f) not in configured and not any(l == lbl for l, _ in dbs):
            dbs.append((lbl, str(f)))
            log.info("homology: auto-including %s from %s", fname, dbdir)
    # Skip DBs whose FASTA isn't actually present (e.g. COG/Pfam listed in config but not yet
    # downloaded) - warn and continue rather than crashing at diamond makedb, so a missing
    # optional database never fails the whole run.
    present = []
    for label, path in dbs:
        if path and Path(path).exists():
            present.append((label, path))
        else:
            log.warning("protein DB '%s' (%s) not found - skipping it in the homology "
                        "filter (get it with download_databases.sh, or fix the path).",
                        label, path)
    return present


def run_protein_filter(ctx: RunContext, tx_fasta: Path,
                       candidates: Set[str]) -> Optional[FilterResult]:
    if not candidates:
        return _empty_result("protein")
    dbs = _protein_databases(ctx)
    if not dbs:
        log.warning("protein-homology filter skipped: set reference.protein_db "
                    "(e.g. uniprot_sprot.fasta) and/or reference.genome_protein "
                    "(your species' proteins) to remove protein-coding transcripts.")
        return None

    backend = ctx.config.tools.get("protein_filter", "diamond")
    outdir = ctx.dirs["protein_filter"]
    fasta = candidate_fasta(ctx, tx_fasta, candidates, "protein")
    evalue = float(ctx.config.get("protein_homology.evalue", 1e-5))
    min_pident = float(ctx.config.get("protein_homology.min_pident", 60))
    min_bit = float(ctx.config.get("protein_homology.min_bitscore", 0))
    threads, logs = ctx.config.threads, ctx.logs

    # Search each configured protein DB (distinct db prefix + hits file per DB so they
    # never collide), and remove a transcript if it matches a protein in ANY of them.
    kept = set(candidates)
    for label, prot_fasta in dbs:
        out_tsv = outdir / f"{backend}_hits_{label}.tsv"
        if backend == "diamond":
            db = outdir / f"diamond_db_{label}"
            if not Path(f"{db}.dmnd").exists():
                get_wrapper("makedb", "diamond")(threads=threads, log_dir=logs).run(
                    protein_fasta=prot_fasta, db_out=db)
            get_wrapper("protein_filter", "diamond")(threads=threads, log_dir=logs).run(
                query_fasta=fasta, db=db, out_tsv=out_tsv, evalue=evalue)
        else:                               # blastx
            db = outdir / f"blast_db_{label}"
            if not Path(f"{db}.pin").exists():
                get_wrapper("makedb", "blast")(threads=threads, log_dir=logs).run(
                    fasta=prot_fasta, db_out=db, dbtype="prot")
            get_wrapper("protein_filter", "blastx")(threads=threads, log_dir=logs).run(
                query_fasta=fasta, db=db, out_tsv=out_tsv, evalue=evalue)
        r = filter_protein_hits(out_tsv, candidates, evalue=evalue,
                                min_pident=min_pident, min_bitscore=min_bit)
        log.info("protein[%s]: %d/%d transcript(s) matched a protein",
                 label, len(candidates) - len(r.kept), len(candidates))
        kept &= r.kept                      # coding if it hits ANY database

    return FilterResult("protein", kept, candidates - kept)


# --------------------------------------------------------------------------- #
# 3. Housekeeping ncRNA removal  (the Plant-LncPipe gap)
# --------------------------------------------------------------------------- #
def run_housekeeping(ctx: RunContext, tx_fasta: Path,
                     candidates: Set[str]) -> Optional[FilterResult]:
    if not candidates:
        return _empty_result("housekeeping")
    backend = ctx.config.tools.get("housekeeping", "infernal")
    if backend in (None, "none"):
        log.warning("housekeeping ncRNA removal is disabled - rRNA/tRNA/snoRNA/miRNA "
                    "will remain in the final set.")
        return None

    hk = ctx.config.get("housekeeping", {}) or {}
    cfg = HousekeepingConfig(
        evalue=float(hk.get("infernal_evalue", 1e-5)),
        remove_biotypes=list(hk.get("remove_biotypes", HousekeepingConfig().remove_biotypes)),
        remove_all_hits=bool(hk.get("remove_all_hits", False)),
        keep_families=list(hk.get("keep_families", [])),
    )
    outdir = ctx.dirs["rfam_filter"]
    fasta = candidate_fasta(ctx, tx_fasta, candidates, "rfam")
    threads, logs = ctx.config.threads, ctx.logs

    if backend == "infernal":
        rfam_cm = ctx.config.get("reference.rfam_cm")
        if not rfam_cm:
            log.warning("housekeeping filter skipped: set reference.rfam_cm to an "
                        "Rfam.cm file (ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/). "
                        "Without it, structural ncRNAs stay in the final set.")
            return None
        rfam_cm = Path(rfam_cm)
        if not Path(f"{rfam_cm}.i1m").exists():
            log.info("Rfam.cm is not indexed - running cmpress (one-off)")
            get_wrapper("cmpress", "infernal")(threads=threads, log_dir=logs).run(
                rfam_cm=rfam_cm)
        tblout = outdir / "rfam_cmscan.tblout"
        get_wrapper("housekeeping", "infernal")(threads=threads, log_dir=logs).run(
            fasta=fasta, rfam_cm=rfam_cm, tblout=tblout, evalue=cfg.evalue,
            use_ga=bool(hk.get("use_rfam_gathering_threshold", True)),
            clanin=hk.get("clanin"))
        hits = parse_cmscan_tblout(tblout, cfg.evalue)
    else:                                    # blastn against an ncRNA FASTA
        ncrna_fasta = ctx.config.get("reference.ncrna_fasta")
        if not ncrna_fasta:
            log.warning("housekeeping filter skipped: tools.housekeeping=blastn needs "
                        "reference.ncrna_fasta (e.g. Rfam.fa or a plant ncRNA set).")
            return None
        db = outdir / "ncrna_db"
        if not Path(f"{db}.nin").exists():
            get_wrapper("makedb", "blast")(threads=threads, log_dir=logs).run(
                fasta=ncrna_fasta, db_out=db, dbtype="nucl")
        out_tsv = outdir / "ncrna_hits.tsv"
        get_wrapper("housekeeping", "blastn")(threads=threads, log_dir=logs).run(
            query_fasta=fasta, db=db, out_tsv=out_tsv, evalue=cfg.evalue)
        hits = parse_blast_hits(out_tsv, cfg.evalue)

    return filter_housekeeping(hits, candidates, cfg)
