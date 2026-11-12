"""DIAMOND wrapper - fast protein homology search (blastx mode).

Transcripts with a significant hit against a protein DB (e.g. UniProt/SwissProt)
are treated as protein-coding and removed from the lncRNA candidate set.
"""
from __future__ import annotations
from pathlib import Path
from typing import List
from lnckit.wrappers.base import BaseWrapper
from lnckit.wrappers.registry import register


@register("makedb", "diamond")
class DiamondMakedbWrapper(BaseWrapper):
    requires = ("diamond",)
    tool_name = "diamond-makedb"

    def build_command(self, *, protein_fasta: Path, db_out: Path, **_) -> List[str]:
        return ["diamond", "makedb", "--in", str(protein_fasta), "-d", str(db_out)]

    def expected_outputs(self, *, db_out: Path, **_) -> List[Path]:
        return [Path(f"{db_out}.dmnd")]


@register("protein_filter", "diamond")
class DiamondBlastxWrapper(BaseWrapper):
    requires = ("diamond",)
    tool_name = "diamond-blastx"
    allow_empty_outputs = True      # zero protein hits is the expected lncRNA result

    def build_command(self, *, query_fasta: Path, db: Path, out_tsv: Path,
                      evalue: float = 1e-5, **_) -> List[str]:
        return ["diamond", "blastx", "-d", str(db), "-q", str(query_fasta),
                "-o", str(out_tsv), "-e", evalue, "-p", self.threads,
                "--outfmt", "6", "qseqid", "sseqid", "pident", "evalue", "bitscore"]

    def expected_outputs(self, *, out_tsv: Path, **_) -> List[Path]:
        return [Path(out_tsv)]
