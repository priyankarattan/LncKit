"""Declarative description of the pipeline DAG.

Kept as data (not hard-coded calls) so the engine can print the plan, support
`--from`/`--until`, skip completed steps on resume, and let plugins insert steps.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Callable, List


@dataclass
class Step:
    name: str
    description: str
    stage: str                       # config key that selects the backend
    produces: List[str] = field(default_factory=list)

# Canonical ordered pipeline (accessions -> upstream -> identification -> report)
PIPELINE: List[Step] = [
    Step("fetch",          "Download SRA accessions -> FASTQ",        "sra"),
    Step("qc",             "Quality control of raw reads",            "qc"),
    Step("trim",           "Adapter/quality trimming",                "trim"),
    Step("merge",          "Pool trimmed replicates by group (cat)",  "-"),
    Step("index",          "Genome indexing",                         "aligner"),
    Step("align",          "Read alignment + BAM sorting",            "aligner"),
    Step("assemble",       "Transcript assembly + merge + FASTA",     "assembler"),
    Step("compare",        "Compare to reference; extract novel",     "compare"),
    Step("filter_length",  "Length filter (>=200 nt)",                "-"),
    Step("filter_exon",    "Exon-number filter",                      "-"),
    Step("filter_expr",    "Expression (TPM) filter",                 "-"),
    Step("filter_orf",     "ORF filter (remove long-ORF transcripts)","orf"),
    Step("coding",         "Coding-potential prediction",             "coding_potential"),
    Step("protein",        "Protein-homology filter",                 "protein_filter"),
    Step("housekeeping",   "Housekeeping ncRNA removal (Rfam)",       "housekeeping"),
    Step("identify",       "Write the high-confidence lncRNA set",    "-"),
    Step("classify",       "Genomic-location classification",         "classifier"),
    Step("report",         "HTML/PDF report + metadata",              "-"),
]

STEP_NAMES = [s.name for s in PIPELINE]
