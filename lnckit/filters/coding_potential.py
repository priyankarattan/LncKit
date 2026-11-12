"""Coding-potential filter - keep only transcripts predicted non-coding.

Backend-agnostic by design: each supported tool has a small parser returning a
per-transcript score or label, so `tools.coding_potential` is swappable in YAML.

Organism-specific models matter here. CPAT and PLEK are trained per organism and
mis-classify at a high rate if the model does not match the species; supply the CPAT
hexamer/logit model for YOUR organism via `reference.cpat_hexamer` / `reference.cpat_logit`
and set `coding_potential.cpat.coding_cutoff` to that model's cutoff (e.g. ~0.36 human,
0.46 plant). CPC2 is model-free and works for any organism as-is; CNCI picks its model
(`-m pl`/`ve`) from the `organism` setting.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Optional, Set

from lnckit.core.logging_utils import get_logger
from lnckit.filters.base import FilterResult

log = get_logger("lnckit.filters.coding")

# CPAT v3 reports one row per ORF, suffixing the transcript id with "_ORF_1", "_ORF_2"...
_ORF_SUFFIX = re.compile(r"_ORF_\d+$")


def _strip_orf_suffix(seq_id: str) -> str:
    return _ORF_SUFFIX.sub("", seq_id)


def _find_col(header, *names) -> Optional[int]:
    lower = [h.strip().lower() for h in header]
    for n in names:
        if n.lower() in lower:
            return lower.index(n.lower())
    return None


# --------------------------------------------------------------------------- #
# CPAT (default backend)
# --------------------------------------------------------------------------- #
def _read_cpat_table(out_prefix: Path) -> Dict[str, Dict[str, float]]:
    """transcript_id -> {'prob': coding_prob, 'orf': ORF size}, max prob across ORFs.

    Handles CPAT v3 (`<prefix>.ORF_prob.best.tsv` + `<prefix>.no_ORF.txt`) and CPAT
    1.2.4's single table. CPAT 1.2.4 writes the sequence id as an UNNAMED first column,
    so its header has one fewer field than the data rows; we detect that shift and offset
    every column index by it - otherwise `coding_prob` (and `ORF_size`) would be read from
    the wrong column. Transcripts CPAT reports with NO ORF are scored prob 0 / orf 0.
    """
    prefix = Path(out_prefix)
    out: Dict[str, Dict[str, float]] = {}

    candidates = [
        prefix.parent / f"{prefix.name}.ORF_prob.best.tsv",   # CPAT v3
        prefix.parent / f"{prefix.name}.ORF_prob.tsv",
        prefix if prefix.is_file() else None,                 # CPAT 1.2.4 single table
        prefix.parent / f"{prefix.name}.dat",
    ]
    table = next((p for p in candidates if p and p.exists()), None)
    if table is None:
        raise FileNotFoundError(
            f"No CPAT result table found for prefix '{prefix}'. Expected "
            f"{prefix.name}.ORF_prob.best.tsv (CPAT v3) or {prefix.name} (CPAT 1.2.4)."
        )

    with open(table) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        first = fh.readline()
        prob_h = _find_col(header, "Coding_prob", "coding_prob", "coding_probability")
        if prob_h is None:
            prob_h = len(header) - 1          # CPAT always puts the probability last
        orf_h = _find_col(header, "ORF_size", "ORF")   # 1.2.4: ORF_size; v3: ORF (nt)
        id_h = _find_col(header, "ID", "seq_ID")

        def emit(line: str, shift: int):
            f = line.rstrip("\n").split("\t")
            pi = prob_h + shift
            if len(f) <= pi:
                return
            try:
                prob = float(f[pi])
            except ValueError:
                return
            tid = _strip_orf_suffix(f[(id_h + shift) if id_h is not None else 0])
            orf = None
            if orf_h is not None and len(f) > orf_h + shift:
                try:
                    orf = float(f[orf_h + shift])
                except ValueError:
                    orf = None
            rec = out.setdefault(tid, {"prob": 0.0, "orf": 0.0})
            if prob >= rec["prob"]:           # keep the ORF of the highest-prob row
                rec["prob"] = prob
                if orf is not None:
                    rec["orf"] = orf

        if first:
            # shift = 1 when the id column is unnamed (data has one more field than header)
            shift = 1 if len(first.rstrip("\n").split("\t")) == len(header) + 1 else 0
            emit(first, shift)
            for line in fh:
                emit(line, shift)

    no_orf = prefix.parent / f"{prefix.name}.no_ORF.txt"
    if no_orf.exists():
        for line in no_orf.read_text().splitlines():
            parts = line.strip().split()
            if parts:
                out.setdefault(_strip_orf_suffix(parts[0]), {"prob": 0.0, "orf": 0.0})
    return out


def parse_cpat(out_prefix: Path) -> Dict[str, float]:
    """transcript_id -> coding probability (kept for callers that only need the prob)."""
    return {t: v["prob"] for t, v in _read_cpat_table(out_prefix).items()}


def filter_cpat(out_prefix: Path, candidates: Set[str],
                coding_cutoff: float = 0.46,
                max_orf_size: Optional[float] = None) -> FilterResult:
    """Keep transcripts with coding_prob < cutoff (non-coding). When max_orf_size is set
    (e.g. the ricebean paper's ORF_size < 100), ALSO require a short ORF - a transcript is
    kept only if BOTH hold. Transcripts CPAT never scored carry no evidence -> kept."""
    table = _read_cpat_table(out_prefix)

    def keep(t: str) -> bool:
        rec = table.get(t)
        if rec is None:
            return True                        # unscored -> no evidence against it
        if rec["prob"] >= coding_cutoff:
            return False
        if max_orf_size is not None and rec["orf"] > max_orf_size:
            return False                       # ORF longer than the cap -> putative coding, drop
        return True

    kept = {t for t in candidates if keep(t)}
    extra = f" and ORF size <= {max_orf_size:g}" if max_orf_size is not None else ""
    log.info("CPAT: %d/%d candidates non-coding (coding_prob < %.2f%s)",
             len(kept), len(candidates), coding_cutoff, extra)
    return FilterResult("coding_potential", kept, candidates - kept)


# --------------------------------------------------------------------------- #
# Label-based backends: CPC2 / PLEK
# --------------------------------------------------------------------------- #
def parse_cpc2(table: Path) -> Dict[str, str]:
    """transcript_id -> 'coding' | 'noncoding' from CPC2's result table."""
    labels: Dict[str, str] = {}
    with open(table) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        lab_i = _find_col(header, "label")
        lab_i = -1 if lab_i is None else lab_i
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            labels[f[0].split()[0]] = f[lab_i].strip().lower()
    return labels


def parse_plek(out_file: Path) -> Dict[str, str]:
    """PLEK writes: '<Coding|Non-coding> <tab> <score> <tab> >seq header'."""
    labels: Dict[str, str] = {}
    with open(out_file) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            tid = f[2].lstrip(">").split()[0]
            labels[tid] = "noncoding" if "non" in f[0].lower() else "coding"
    return labels


def parse_pfam(domtblout: Path, candidates: Set[str],
               evalue: float = 1e-5) -> Dict[str, str]:
    """Label transcripts coding if their translated ORF hit a Pfam domain.

    Reads an HMMER --domtblout table. The getorf peptide id ('<transcript>_<n>') is the
    'query' (col 4) for hmmscan but the 'target' (col 1) for hmmsearch - the two tools
    swap those columns. We auto-detect which column holds the peptide by matching it back
    to `candidates`, so the same parser works for EITHER tool. Column 7 is the
    full-sequence E-value in both. A transcript is 'coding' if ANY of its ORF peptides has
    a domain hit with E <= evalue; getorf appends '_<number>' per ORF, which we strip to
    recover the transcript id. Any candidate without a qualifying hit stays 'noncoding'.
    """
    labels: Dict[str, str] = {t: "noncoding" for t in candidates}
    try:
        fh = open(domtblout)
    except OSError:
        return labels

    def _to_tid(name: str) -> Optional[str]:
        """Return the candidate transcript id this domtbl name belongs to, or None."""
        base = re.sub(r"_\d+$", "", name)      # drop getorf's ORF suffix
        if base in labels:
            return base
        if name in labels:                     # (defensive: no suffix was added)
            return name
        return None

    with fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 7:
                continue
            try:
                ev = float(f[6])               # full-sequence E-value (same col both tools)
            except ValueError:
                continue
            if ev > evalue:
                continue
            # The peptide is col 4 (hmmscan) or col 1 (hmmsearch); take whichever maps to
            # a candidate. The other column is the Pfam profile name/accession.
            tid = _to_tid(f[3]) or _to_tid(f[0])
            if tid is not None:
                labels[tid] = "coding"
    return labels


def parse_rnaplonc(cds_file: Path, weka_out: Path,
                   candidates: Set[str]) -> Dict[str, str]:
    """Label transcripts from RNAplonc's Weka REPTree output.

    RNAplonc emits one CDS row per transcript (column 0 = transcript id) and feeds them,
    in that same order, through feature_extraction.pl -> Weka. Weka `-p 0` then prints one
    prediction per row, in order, as e.g. `  1   1:?   1:lncRNA   0.863`. We therefore
    zip the CDS ids (row order) with the Weka predictions (row order): '1:lncRNA' ->
    non-coding (kept), '2:mRNA' -> coding. Unmatched candidates default to non-coding.
    """
    labels: Dict[str, str] = {t: "noncoding" for t in candidates}
    # CDS ids in file order (column 0)
    ids = []
    try:
        with open(cds_file) as fh:
            for line in fh:
                if line.strip():
                    ids.append(line.split("\t")[0].strip())
    except OSError:
        return labels
    # Weka predictions in order (start after the 'inst# ... predicted ...' header)
    preds = []
    started = False
    try:
        with open(weka_out) as fh:
            for line in fh:
                if not started:
                    if "inst#" in line and "predicted" in line:
                        started = True
                    continue
                if not line.strip():
                    continue
                # token like '2:mRNA' or '1:lncRNA' is the predicted class
                m = re.search(r"\d+:(lncRNA|mRNA)", line)
                if m:
                    preds.append(m.group(1))
    except OSError:
        return labels
    for tid, pred in zip(ids, preds):
        if tid in labels:
            labels[tid] = "coding" if pred == "mRNA" else "noncoding"
    return labels


def parse_lncfinder(tsv: Path) -> Dict[str, str]:
    """transcript_id -> 'coding' | 'noncoding' from the LncFinder driver's TSV.

    LncFinder's Pred column holds "Coding" / "NonCoding"; "non" identifies the latter.
    """
    labels: Dict[str, str] = {}
    with open(tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        pred_i = _find_col(header, "Pred", "pred")
        pred_i = 1 if pred_i is None else pred_i
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= pred_i:
                continue
            tid = f[0].strip().strip('"').split()[0]
            labels[tid] = "noncoding" if "non" in f[pred_i].lower() else "coding"
    return labels


def parse_cnci(index_file: Path) -> Dict[str, str]:
    """transcript_id -> 'coding' | 'noncoding' from CNCI's `CNCI.index` table.

    CNCI writes a header row then one line per transcript; the second column
    ("index") holds 'coding' or 'noncoding'. The transcript id in the first column
    keeps CNCI's leading '>' and any trailing description, so it is normalised here.
    """
    labels: Dict[str, str] = {}
    with open(index_file) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        lab_i = _find_col(header, "index", "label")
        lab_i = 1 if lab_i is None else lab_i
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= lab_i:
                continue
            tid = f[0].lstrip(">").split()[0]
            labels[tid] = "noncoding" if "non" in f[lab_i].lower() else "coding"
    return labels


def filter_by_label(labels: Dict[str, str], candidates: Set[str],
                    tool: str = "cpc2") -> FilterResult:
    """Keep candidates the tool labelled non-coding (unscored transcripts are kept)."""
    kept = {t for t in candidates if labels.get(t, "noncoding") == "noncoding"}
    log.info("%s: %d/%d candidates labelled non-coding", tool, len(kept), len(candidates))
    return FilterResult("coding_potential", kept, candidates - kept)


def filter_cpc2(cpc2_table: Path, candidates: Set[str]) -> FilterResult:
    return filter_by_label(parse_cpc2(cpc2_table), candidates, tool="cpc2")
