"""STAR wrapper - fast splice-aware aligner (higher RAM than HISAT2)."""
from __future__ import annotations
from pathlib import Path
from typing import List, Optional
from lnckit.wrappers.base import BaseWrapper
from lnckit.wrappers.registry import register


@register("index", "star")
class StarIndexWrapper(BaseWrapper):
    requires = ("star",)
    tool_name = "STAR-index"

    def build_command(self, *, genome: Path, index_dir: Path,
                      annotation: Optional[Path] = None, sjdb_overhang: int = 100, **_) -> List[str]:
        cmd = ["STAR", "--runMode", "genomeGenerate", "--runThreadN", self.threads,
               "--genomeDir", str(index_dir), "--genomeFastaFiles", str(genome)]
        if annotation:
            cmd += ["--sjdbGTFfile", str(annotation), "--sjdbOverhang", sjdb_overhang]
        return cmd

    def expected_outputs(self, *, index_dir: Path, **_) -> List[Path]:
        return [Path(index_dir) / "SAindex"]


@register("aligner", "star")
class StarAlignWrapper(BaseWrapper):
    requires = ("star",)
    tool_name = "STAR"

    def build_command(self, *, index_dir: Path, out_prefix: Path,
                      fq1: Path, fq2: Optional[Path] = None, **_) -> List[str]:
        reads = [str(fq1)] + ([str(fq2)] if fq2 else [])
        cmd = ["STAR", "--runThreadN", self.threads, "--genomeDir", str(index_dir),
               "--readFilesIn", *reads]
        # STAR does NOT auto-decompress: gzipped reads need an explicit reader, or STAR
        # reads the gzip bytes as sequence and finds ~0 alignments.
        if str(fq1).endswith(".gz"):
            cmd += ["--readFilesCommand", "zcat"]
        # --outSAMstrandField intronMotif adds the XS:A: strand tag to spliced alignments
        # (from the splice-site motif). StringTie/Cufflinks need it for STAR output, and it
        # is what the strandedness auto-detect reads; harmless for already-stranded data.
        cmd += ["--outSAMtype", "BAM", "SortedByCoordinate",
                "--outSAMstrandField", "intronMotif",
                "--outFileNamePrefix", str(out_prefix)]
        return cmd

    def expected_outputs(self, *, out_prefix: Path, **_) -> List[Path]:
        return [Path(f"{out_prefix}Aligned.sortedByCoord.out.bam")]
