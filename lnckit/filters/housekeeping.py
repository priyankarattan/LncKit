"""Housekeeping-ncRNA filter - the step Plant-LncPipe leaves out.

Coding-potential models and protein-homology search both answer one question:
"could this transcript encode a protein?" A tRNA, an rRNA fragment, a snoRNA or a
miRNA precursor answers "no" to that question just as emphatically as a real
lncRNA does, so every upstream filter *passes them straight through*. Without this
stage the final set is "non-coding RNA", not "long non-coding RNA".

This module removes any candidate with a significant hit to a structural /
housekeeping ncRNA family, parsed from either:
  * Infernal `cmscan --tblout` (covariance models; both --fmt 1 and --fmt 2), or
  * BLASTN outfmt-6 against an ncRNA FASTA.

Family assignment
-----------------
Rfam family *names* are not self-describing (`U1` is a snRNA, `SNORD13` a snoRNA,
`MIR166` a miRNA), so each biotype is matched against the family name **and** the
target description. `remove_all_hits` implements the stricter policy of dropping
every Rfam hit; because Rfam also curates genuine lncRNA families, anything listed
in `keep_families` is protected from removal.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set

from lnckit.core.logging_utils import get_logger
from lnckit.filters.base import FilterResult

log = get_logger("lnckit.filters.housekeeping")

# Rfam curates genuine lncRNA families (e.g. XIST, HOTAIR, MALAT1, lincRNAs,
# antisense RNAs). A transcript whose Rfam hit is such a family must NEVER be removed
# by the housekeeping step - that would throw away exactly what we are trying to find.
# This guard is applied unconditionally, before any biotype/remove_all_hits rule.
_LNCRNA_PROTECT = re.compile(
    r"lncrna|long non[- ]?coding|antisense|\blinc[-_ ]?rna|\bxist\b|\bhotair\b|\bmalat", re.I)

# Biotype -> pattern matched (case-insensitively) against "family name + description".
BIOTYPE_PATTERNS: Dict[str, str] = {
    "rrna":       r"rrna|ribosomal rna|\blsu\b|\bssu\b|\b5s\b|\b5_8s\b|\b16s\b|\b18s\b|\b23s\b|\b25s\b|\b28s\b",
    "trna":       r"\btrna\b|transfer rna|trna-",
    "snrna":      r"\bsnrna\b|small nuclear|spliceosomal|^u[1-6](atac)?$|^u1[12]$",
    "snorna":     r"snorna|small nucleolar|^snor|^sno[a-z]|scarna|\bh/aca\b|\bc/d\b",
    "mirna":      r"\bmirna\b|microrna|^mir[-_]?\d|^mir\b",
    "srprna":     r"\bsrp\b|signal recognition",
    "rnasep":     r"rnase[ _]?p\b|rnase[ _]?mrp",
    "tmrna":      r"\btmrna\b",
    "ribozyme":   r"ribozyme|hammerhead",
    "riboswitch": r"riboswitch",
    "yrna":       r"\by[ _]?rna\b",
    "telomerase": r"telomerase",
    "vault":      r"vault",
}


@dataclass
class HousekeepingConfig:
    evalue: float = 1e-5
    remove_biotypes: List[str] = field(default_factory=lambda: [
        "rRNA", "tRNA", "snRNA", "snoRNA", "miRNA", "srpRNA", "rnasep", "tmRNA",
    ])
    #: drop EVERY significant Rfam hit, not only the listed biotypes
    remove_all_hits: bool = False
    #: families protected from removal (Rfam curates real lncRNA families too)
    keep_families: List[str] = field(default_factory=list)


def _compiled(biotypes: Iterable[str]) -> List[re.Pattern]:
    pats = []
    for b in biotypes:
        key = b.strip().lower().replace("-", "").replace("_", "")
        pattern = BIOTYPE_PATTERNS.get(key)
        # an unrecognised biotype is matched literally against name/description
        pats.append(re.compile(pattern or re.escape(b.strip()), re.I))
    return pats


def parse_cmscan_tblout(tblout: Path, evalue: float = 1e-5) -> Dict[str, Set[str]]:
    """Return transcript_id -> {"family: description"} for significant Rfam hits.

    Handles both Infernal tblout layouts, detected by whether the first field is a
    hit index:
      --fmt 1 : target(0) accession(1) query(2) ... E-value(15) ... description(17+)
      --fmt 2 : idx(0) target(1) accession(2) query(3) ... E-value(17) olp(19) ... description(26+)
    """
    hits: Dict[str, Set[str]] = {}
    with open(tblout) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            try:
                if f[0].isdigit():                       # --fmt 2
                    # '=' in the overlap column marks a hit already represented by a
                    # better-scoring one in the same Rfam clan.
                    if len(f) > 19 and f[19] == "=":
                        continue
                    family, query, ev, desc_at = f[1], f[3], float(f[17]), 26
                else:                                    # --fmt 1
                    family, query, ev, desc_at = f[0], f[2], float(f[15]), 17
            except (IndexError, ValueError):
                continue
            if ev > evalue:
                continue
            description = " ".join(f[desc_at:]) if len(f) > desc_at else ""
            hits.setdefault(query, set()).add(f"{family}: {description}".strip())
    return hits


def parse_blast_hits(tsv: Path, evalue: float = 1e-5) -> Dict[str, Set[str]]:
    """transcript_id -> {subject_id} from BLASTN outfmt 6 (qseqid sseqid pident evalue bitscore)."""
    hits: Dict[str, Set[str]] = {}
    with open(tsv) as fh:
        for line in fh:
            f = line.rstrip().split("\t")
            if len(f) < 4:
                continue
            try:
                ev = float(f[3])
            except ValueError:
                continue
            if ev <= evalue:
                hits.setdefault(f[0], set()).add(f[1])
    return hits


def filter_housekeeping(hits: Dict[str, Set[str]], candidates: Set[str],
                        cfg: Optional[HousekeepingConfig] = None) -> FilterResult:
    """Remove candidates whose ncRNA hits identify them as housekeeping ncRNA."""
    cfg = cfg or HousekeepingConfig()
    patterns = _compiled(cfg.remove_biotypes)
    protect = _compiled(cfg.keep_families) if cfg.keep_families else []

    removed: Set[str] = set()
    protected_lnc = 0
    for tid in candidates:
        annotations = hits.get(tid)
        if not annotations:
            continue
        blob = " ".join(annotations)
        # Never remove a transcript that hits a real lncRNA family, regardless of policy.
        if _LNCRNA_PROTECT.search(blob):
            protected_lnc += 1
            continue
        if protect and any(p.search(blob) for p in protect):
            continue
        if cfg.remove_all_hits or any(p.search(blob) for p in patterns):
            removed.add(tid)

    if removed:
        log.info("housekeeping filter removed %d structural ncRNA transcript(s) "
                 "(rRNA/tRNA/snRNA/snoRNA/miRNA/...)", len(removed))
    if protected_lnc:
        log.info("housekeeping filter kept %d transcript(s) that hit an Rfam lncRNA "
                 "family (protected from removal)", protected_lnc)
    return FilterResult("housekeeping", candidates - removed, removed)


def filter_infernal(tblout: Path, candidates: Set[str], evalue: float = 1e-5,
                    remove_biotypes: Iterable[str] = ()) -> FilterResult:
    """Convenience entry point: parse a cmscan tblout and filter in one call."""
    cfg = HousekeepingConfig(evalue=evalue)
    if remove_biotypes:
        cfg.remove_biotypes = list(remove_biotypes)
    return filter_housekeeping(parse_cmscan_tblout(tblout, evalue), candidates, cfg)
