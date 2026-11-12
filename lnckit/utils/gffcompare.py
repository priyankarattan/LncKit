"""Parsing of gffcompare output - how novel transcripts are separated from known ones.

gffcompare assigns every assembled transcript a *class code* describing its
relationship to the reference annotation. The ones that matter for lncRNA discovery:

    u   intergenic - no overlap with any reference gene   (lincRNA candidates)
    x   exonic overlap on the *opposite* strand           (antisense candidates)
    i   fully contained within a reference intron         (intronic candidates)
    o   generic exonic overlap on the same strand         (sense-overlapping)
    =   exact match to a reference transcript             -> known, discarded
    j   novel isoform of a reference gene                 -> usually protein-coding, discarded
    c   contained within a reference transcript           -> discarded

Keeping only {u, x, i, o} (configurable via `filters.keep_classes`) removes the
bulk of protein-coding assemblies before any expensive model is run.

The `.tmap` file is one row per *query* transcript, which is exactly the granularity
the filters work at:
    ref_gene_id  ref_id  class_code  qry_gene_id  qry_id  num_exons  FPKM  TPM  cov  len
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set

from lnckit.core.logging_utils import get_logger

log = get_logger("lnckit.gffcompare")


@dataclass
class TmapEntry:
    qry_id: str
    class_code: str
    ref_gene_id: str = "-"
    ref_id: str = "-"
    num_exons: int = 0
    tpm: float = 0.0
    length: int = 0


def find_tmap(compare_prefix: Path) -> Optional[Path]:
    """gffcompare names it '<prefix>.<query-gtf-basename>.tmap'."""
    prefix = Path(compare_prefix)
    hits = sorted(prefix.parent.glob(f"{prefix.name}*.tmap"))
    return hits[0] if hits else None


def parse_tmap(tmap: Path) -> Dict[str, TmapEntry]:
    """Return qry_transcript_id -> TmapEntry."""
    entries: Dict[str, TmapEntry] = {}
    with open(tmap) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}

        def get(fields: List[str], name: str, default: str = "") -> str:
            i = col.get(name, -1)
            return fields[i] if 0 <= i < len(fields) else default

        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            qry = get(f, "qry_id")
            if not qry:
                continue
            entries[qry] = TmapEntry(
                qry_id=qry,
                class_code=get(f, "class_code", "-"),
                ref_gene_id=get(f, "ref_gene_id", "-"),
                ref_id=get(f, "ref_id", "-"),
                num_exons=int(get(f, "num_exons", "0") or 0),
                tpm=float(get(f, "TPM", "0") or 0),
                length=int(get(f, "len", "0") or 0),
            )
    return entries


def novel_transcript_ids(tmap: Path, keep_classes: Iterable[str]) -> Set[str]:
    """Transcript IDs whose class code marks them as novel lncRNA candidates."""
    keep = {c.strip() for c in keep_classes}
    entries = parse_tmap(tmap)
    novel = {e.qry_id for e in entries.values() if e.class_code in keep}
    by_code: Dict[str, int] = {}
    for e in entries.values():
        by_code[e.class_code] = by_code.get(e.class_code, 0) + 1
    log.info("gffcompare class codes: %s",
             ", ".join(f"{k}={v}" for k, v in sorted(by_code.items())))
    log.info("kept %d/%d transcripts with class codes %s",
             len(novel), len(entries), sorted(keep))
    return novel
