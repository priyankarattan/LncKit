"""NCBI BLAST+ wrappers - protein-homology (blastx) and ncRNA (blastn) backends.

BLASTX is the classical alternative to DIAMOND for the protein-homology filter:
far slower, but the reference implementation many reviewers expect. BLASTN against
an Rfam/ncRNA FASTA is the lightweight alternative to Infernal for housekeeping-ncRNA
removal when no covariance-model search is available (sequence-only, so less
sensitive to the structural families Infernal catches).

Both emit outfmt 6 with the same five columns the filters parse:
    qseqid sseqid pident evalue bitscore
"""
from __future__ import annotations
from pathlib import Path
from typing import List

from lnckit.wrappers.base import BaseWrapper
from lnckit.wrappers.registry import register

OUTFMT6 = "6 qseqid sseqid pident evalue bitscore"


@register("makedb", "blast")
class MakeBlastDbWrapper(BaseWrapper):
    requires = ("makeblastdb",)
    tool_name = "makeblastdb"

    def build_command(self, *, fasta: Path, db_out: Path, dbtype: str = "prot", **_) -> List[str]:
        return ["makeblastdb", "-in", str(fasta), "-dbtype", dbtype, "-out", str(db_out)]

    def expected_outputs(self, *, db_out: Path, dbtype: str = "prot", **_) -> List[Path]:
        # BLAST+ >=2.10 writes .p?? / .n?? files; the index file is the stable one.
        suffix = "pin" if dbtype == "prot" else "nin"
        return [Path(f"{db_out}.{suffix}")]


@register("protein_filter", "blastx")
class BlastxWrapper(BaseWrapper):
    requires = ("blastx",)
    tool_name = "blastx"
    allow_empty_outputs = True      # "no protein hit" is the expected result here

    def build_command(self, *, query_fasta: Path, db: Path, out_tsv: Path,
                      evalue: float = 1e-5, **_) -> List[str]:
        return ["blastx", "-query", str(query_fasta), "-db", str(db),
                "-out", str(out_tsv), "-evalue", evalue,
                "-num_threads", self.threads, "-max_target_seqs", 1,
                "-outfmt", OUTFMT6]

    def expected_outputs(self, *, out_tsv: Path, **_) -> List[Path]:
        return [Path(out_tsv)]


@register("housekeeping", "blastn")
class BlastnWrapper(BaseWrapper):
    requires = ("blastn",)
    tool_name = "blastn"
    allow_empty_outputs = True      # no ncRNA hit = nothing to remove

    def build_command(self, *, query_fasta: Path, db: Path, out_tsv: Path,
                      evalue: float = 1e-5, **_) -> List[str]:
        return ["blastn", "-query", str(query_fasta), "-db", str(db),
                "-out", str(out_tsv), "-evalue", evalue,
                "-num_threads", self.threads, "-max_target_seqs", 5,
                "-outfmt", OUTFMT6]

    def expected_outputs(self, *, out_tsv: Path, **_) -> List[Path]:
        return [Path(out_tsv)]
