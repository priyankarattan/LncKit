"""Expression filter - drop transcripts that are not credibly expressed anywhere.

lncRNAs are, as a class, expressed at low levels and in a strongly tissue- or
condition-specific way. Averaging TPM across samples therefore *systematically*
discards real lncRNAs: a transcript at 8 TPM in one of twelve libraries averages
below a 0.5 cutoff and disappears.

The test is therefore applied per sample and then counted: a transcript is kept if
it clears `min_tpm` in at least `min_samples` libraries (default 1). With a single
sample this reduces to the familiar `TPM >= min_tpm`.
"""
from __future__ import annotations
from typing import Dict, List, Mapping, Set

from lnckit.filters.base import FilterResult


def filter_by_expression(tpm: Mapping[str, object], candidates: Set[str],
                         min_tpm: float = 0.5, min_samples: int = 1) -> FilterResult:
    """`tpm` maps transcript_id -> float (one sample) or -> list[float] (per sample)."""
    kept: Set[str] = set()
    for t in candidates:
        values = tpm.get(t, 0.0)
        if isinstance(values, (int, float)):
            values = [float(values)]
        n_expressed = sum(1 for v in values if float(v) >= min_tpm)
        if n_expressed >= max(1, min_samples):
            kept.add(t)
    return FilterResult("expression", kept, candidates - kept)


def merge_sample_tpms(per_sample: Dict[str, Dict[str, float]]) -> Dict[str, List[float]]:
    """{sample: {tid: tpm}} -> {tid: [tpm, ... one per sample]} (missing = 0.0)."""
    samples = list(per_sample)
    all_ids: Set[str] = set()
    for table in per_sample.values():
        all_ids.update(table)
    return {tid: [per_sample[s].get(tid, 0.0) for s in samples] for tid in all_ids}
