"""Classify final lncRNAs by genomic position relative to reference genes.

Categories follow the conventional scheme (intergenic/lincRNA, intronic,
antisense, sense-overlapping, upstream/downstream, bidirectional). Two backends:
FEELnc_classifier (reference-based) or gffcompare class codes.
"""
from __future__ import annotations
import re
from pathlib import Path
from typing import Optional

from lnckit.core.context import RunContext
from lnckit.core.logging_utils import get_logger
from lnckit.utils.gffcompare import find_tmap, parse_tmap
from lnckit.wrappers import get_wrapper

log = get_logger("lnckit.classify")

_TID_RE = re.compile(r'transcript_id "([^"]+)"')

# gffcompare class code -> lncRNA positional category (see utils/gffcompare docstring).
_GFFCOMPARE_CLASS = {
    "u": "intergenic",         # no overlap with any reference gene (lincRNA)
    "x": "antisense",          # exonic overlap on the opposite strand
    "s": "antisense",          # intron overlap on the opposite strand
    "i": "intronic",           # fully within a reference intron
    "o": "sense_overlapping",  # generic same-strand exonic overlap
    "e": "sense_overlapping",  # single-exon transfrag overlapping a reference exon
    "p": "downstream",         # polymerase run-on: near a gene, not overlapping
    "=": "reference_match",    # matches/contained in a reference transcript (rare in a
    "c": "reference_match",    # final lncRNA set, but reported rather than silently dropped)
    "j": "reference_match",
    "k": "reference_match",
    "m": "reference_match",
    "n": "reference_match",
    "y": "reference_match",
}


def classify(ctx: RunContext, lncrna_gtf: Path,
             compare_prefix: Optional[Path] = None) -> Path:
    backend = ctx.config.tools.get("classifier", "feelnc")
    ann = ctx.config.get("reference.annotation_gtf")
    out = ctx.dirs["classification"] / "lncRNA_classes.txt"
    if backend == "feelnc":
        w = get_wrapper("classifier", "feelnc")(threads=ctx.config.threads, log_dir=ctx.logs)
        w.run(stdout_path=out, lncrna_gtf=lncrna_gtf, annotation=ann)
    else:
        _classify_by_gffcompare(ctx, lncrna_gtf, compare_prefix, out)
    return out


def _transcript_ids(gtf: Path):
    """Unique transcript ids in GTF order (keeps the output table stable/reproducible)."""
    ids, seen = [], set()
    for line in Path(gtf).read_text().splitlines():
        if line.startswith("#"):
            continue
        m = _TID_RE.search(line)
        if m and m.group(1) not in seen:
            seen.add(m.group(1))
            ids.append(m.group(1))
    return ids


def _classify_by_gffcompare(ctx: RunContext, lncrna_gtf: Path,
                            compare_prefix: Optional[Path], out: Path) -> None:
    """Assign each final lncRNA a positional class from the gffcompare .tmap class codes.

    The compare stage already ran gffcompare against the reference annotation and wrote a
    `.tmap` (one row per assembled transcript). The final lncRNA ids are a subset of those
    query ids, so we look each one up and map its class code to a category. Writes a
    `transcript_id<TAB>class_code<TAB>class` table - the same granularity FEELnc produces.
    """
    prefix = compare_prefix or (ctx.tmp / "gffcompare" / "gffcmp")
    tmap = find_tmap(Path(prefix)) if prefix else None
    if not tmap:
        log.warning("gffcompare classification: no .tmap found beside %s - the compare step "
                    "may not have run. Writing an empty class table; use classifier: feelnc "
                    "for reference-based classification.", prefix)
        out.write_text("transcript_id\tclass_code\tclass\n")
        return

    entries = parse_tmap(tmap)
    counts: dict = {}
    with open(out, "w") as fh:
        fh.write("transcript_id\tclass_code\tclass\n")
        for tid in _transcript_ids(lncrna_gtf):
            e = entries.get(tid)
            code = e.class_code if e else "-"
            category = _GFFCOMPARE_CLASS.get(code, "other")
            counts[category] = counts.get(category, 0) + 1
            fh.write(f"{tid}\t{code}\t{category}\n")
    log.info("gffcompare classification (%s): %s", tmap.name,
             ", ".join(f"{k}={v}" for k, v in sorted(counts.items())) or "no transcripts")
