"""ORF filter - remove transcripts with a long open reading frame.

This is a core lncRNA discriminator. A transcript is flagged as *putative coding*
(and removed from the lncRNA set) if EITHER:
  * its longest ORF exceeds `max_orf_length_nt` (or the aa-equivalent), OR
  * the ORF covers more than `max_orf_coverage` of the transcript length.

Two evidence sources are supported:
  1. External ORF caller output (getorf) parsed for the per-transcript max ORF.
  2. A pure-Python longest-ORF fallback (utils.fasta.longest_orf_length) when no
     external caller is available - keeps the filter usable in CI and tests.
"""
from __future__ import annotations
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Set

from lnckit.core.logging_utils import get_logger
from lnckit.filters.base import FilterResult
from lnckit.utils.fasta import read_fasta, longest_orf_length

log = get_logger("lnckit.filters.orf")

# getorf headers look like: ">TRANSCRIPT_1 [12 - 350] ..." (coords give ORF length)
_GETORF_HDR = re.compile(r">(\S+?)_\d+\s+\[(\d+)\s*-\s*(\d+)\]")


@dataclass
class OrfConfig:
    enabled: bool = True
    max_orf_length_nt: int = 300
    max_orf_length_aa: int = 100
    max_orf_coverage: float = 0.60
    min_transcript_for_orf: int = 200
    strand: str = "both"

    @property
    def effective_max_nt(self) -> int:
        # Use the stricter of the nt limit and the aa-derived nt limit.
        return min(self.max_orf_length_nt, self.max_orf_length_aa * 3)


def parse_getorf_max_orf(getorf_fasta: Path) -> Dict[str, int]:
    """Return transcript_id -> longest ORF length (nt) from getorf output."""
    best: Dict[str, int] = {}
    for name_line, _seq in read_fasta(getorf_fasta):
        # read_fasta strips at first whitespace, so re-read raw headers instead:
        pass
    # getorf headers carry coordinates only in the raw '>' line, so scan the file.
    with open(getorf_fasta) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            m = _GETORF_HDR.search(line)
            if not m:
                continue
            tid, a, b = m.group(1), int(m.group(2)), int(m.group(3))
            length = abs(b - a) + 1
            if length > best.get(tid, 0):
                best[tid] = length
    return best


def compute_orf_lengths_python(fasta: Path, candidates: Set[str], strand: str) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    for name, seq in read_fasta(fasta):
        if name in candidates:
            lengths[name] = longest_orf_length(seq, strand=strand)
    return lengths


def filter_by_orf(
    fasta: Path,
    candidates: Set[str],
    cfg: OrfConfig,
    getorf_output: Optional[Path] = None,
) -> FilterResult:
    """Keep transcripts whose longest ORF is short enough to be non-coding."""
    if not cfg.enabled:
        return FilterResult("orf", set(candidates))

    # transcript lengths (for coverage test)
    tx_len = {n: len(s) for n, s in read_fasta(fasta) if n in candidates}

    if getorf_output and Path(getorf_output).exists():
        orf_len = parse_getorf_max_orf(getorf_output)
        source = "getorf"
    else:
        orf_len = compute_orf_lengths_python(fasta, candidates, cfg.strand)
        source = "python-fallback"
    log.info("ORF filter using %s; max-ORF ceiling=%d nt, coverage<=%.2f",
             source, cfg.effective_max_nt, cfg.max_orf_coverage)

    kept: Set[str] = set()
    for t in candidates:
        L = tx_len.get(t, 0)
        if L < cfg.min_transcript_for_orf:
            # too short to judge by ORF; length filter handles these
            kept.add(t)
            continue
        orf = orf_len.get(t, 0)
        coverage = (orf / L) if L else 0.0
        is_coding = (orf > cfg.effective_max_nt) or (coverage > cfg.max_orf_coverage)
        if not is_coding:
            kept.add(t)
    return FilterResult("orf", kept, candidates - kept)
