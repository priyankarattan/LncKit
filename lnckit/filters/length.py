"""Length filter - enforce the canonical lncRNA >=200 nt floor."""
from __future__ import annotations
from pathlib import Path
from typing import Set
from lnckit.filters.base import FilterResult
from lnckit.utils.fasta import fasta_lengths


def filter_by_length(fasta: Path, candidates: Set[str], min_length: int = 200) -> FilterResult:
    lengths = fasta_lengths(fasta)
    kept = {t for t in candidates if lengths.get(t, 0) >= min_length}
    return FilterResult("length", kept, candidates - kept)
