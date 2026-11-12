"""Protein-homology filter - remove transcripts with a significant protein hit.

A transcript that aligns (DIAMOND/BLASTX) to a known protein is almost certainly
protein-coding and is removed from the lncRNA set.

Thresholds follow the Plant-LncPipe Swissprot recipe: a hit counts as coding when
its **percent identity > 60** and its **e-value < 1e-5**. Percent identity - not
bitscore - is the discriminating criterion there, so it is the primary gate here;
`min_bitscore` is available as an optional extra constraint (0 = off).

Expected input: BLAST/DIAMOND outfmt 6 with columns
    qseqid sseqid pident evalue bitscore
(the 5-column subset lnckit's wrappers emit).
"""
from __future__ import annotations
from pathlib import Path
from typing import Set

from lnckit.core.logging_utils import get_logger
from lnckit.filters.base import FilterResult

log = get_logger("lnckit.filters.protein")


def filter_protein_hits(diamond_tsv: Path, candidates: Set[str],
                        evalue: float = 1e-5, min_pident: float = 60.0,
                        min_bitscore: float = 0.0) -> FilterResult:
    coding: Set[str] = set()
    with open(diamond_tsv) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            qid = parts[0]
            try:
                pident = float(parts[2])
                e = float(parts[3])
                bit = float(parts[4]) if len(parts) > 4 else 0.0
            except ValueError:
                continue
            if e <= evalue and pident >= min_pident and bit >= min_bitscore:
                coding.add(qid)
    kept = {t for t in candidates if t not in coding}
    log.info("protein-homology: %d/%d candidates have no significant protein hit "
             "(pident >= %.0f, e <= %.0e)", len(kept), len(candidates), min_pident, evalue)
    return FilterResult("protein_homology", kept, candidates & coding)
