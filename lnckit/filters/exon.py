"""Exon-number filter - optionally require multi-exon transcripts."""
from __future__ import annotations
from pathlib import Path
from typing import Set
from lnckit.filters.base import FilterResult
from lnckit.utils.gtf import exon_counts


def filter_by_exons(gtf: Path, candidates: Set[str], min_exons: int = 1) -> FilterResult:
    counts = exon_counts(gtf)
    kept = {t for t in candidates if counts.get(t, 0) >= min_exons}
    return FilterResult("exon", kept, candidates - kept)
