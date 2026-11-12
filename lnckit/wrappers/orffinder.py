"""ORF-detection wrappers - central to lncRNA identification.

Biological rationale
--------------------
Protein-coding transcripts contain a long, conserved open reading frame (ORF),
whereas lncRNAs characteristically lack one (typically shorter than ~100 codons
/ 300 nt, and covering a small fraction of the transcript). Computing the longest
ORF per transcript therefore yields a strong, tool-agnostic discriminative
feature that complements probabilistic coding-potential models (CPAT/CPC2).

Three interchangeable backends are provided; all emit a per-transcript longest-ORF
length that filters/orf.py consumes:
  * getorf       (EMBOSS)        - fast, ubiquitous, easy conda install
  * transdecoder (TransDecoder)  - ORF + Pfam-aware, common in assembly pipelines
  * orffinder    (NCBI ORFfinder)- exact NCBI ORF calling
"""
from __future__ import annotations
from pathlib import Path
from typing import List
from lnckit.wrappers.base import BaseWrapper
from lnckit.wrappers.registry import register


@register("orf", "getorf")
class GetorfWrapper(BaseWrapper):
    requires = ("getorf",)
    tool_name = "getorf"

    def build_command(self, *, fasta: Path, orf_out: Path,
                      min_nt: int = 30, genetic_code: int = 1,
                      strand: str = "both", find: int = 0, **_) -> List[str]:
        # -find 0 => nucleotide ORFs between STOP codons (for the ORF-length filter).
        # -find 1 => translated PROTEIN of those ORFs (what the Pfam scan searches).
        # getorf uses -minsize in nucleotides in both cases.
        cmd = ["getorf", "-sequence", str(fasta), "-outseq", str(orf_out),
               "-find", int(find), "-minsize", int(min_nt), "-table", int(genetic_code)]
        if strand == "forward":
            cmd += ["-reverse", "N"]
        return cmd

    def expected_outputs(self, *, orf_out: Path, **_) -> List[Path]:
        return [Path(orf_out)]


@register("orf", "transdecoder")
class TransDecoderWrapper(BaseWrapper):
    requires = ("transdecoder",)
    tool_name = "TransDecoder.LongOrfs"

    def build_command(self, *, fasta: Path, out_dir: Path, min_aa: int = 100, **_) -> List[str]:
        return ["TransDecoder.LongOrfs", "-t", str(fasta), "-m", int(min_aa),
                "--output_dir", str(out_dir)]

    def expected_outputs(self, *, out_dir: Path, **_) -> List[Path]:
        return [Path(out_dir) / "longest_orfs.pep"]


@register("orf", "orffinder")
class OrfFinderWrapper(BaseWrapper):
    requires = ("orffinder",)
    tool_name = "ORFfinder"

    def build_command(self, *, fasta: Path, orf_out: Path,
                      min_nt: int = 30, genetic_code: int = 1, strand: str = "both", **_) -> List[str]:
        strand_flag = {"both": "both", "forward": "plus"}.get(strand, "both")
        return ["ORFfinder", "-in", str(fasta), "-out", str(orf_out),
                "-ml", int(min_nt), "-g", int(genetic_code), "-s", "0",
                "-strand", strand_flag, "-outfmt", "1"]

    def expected_outputs(self, *, orf_out: Path, **_) -> List[Path]:
        return [Path(orf_out)]
