"""Bowtie2 and BWA wrappers - unspliced aligners.

They are supported because users ask for them, but they are **not** valid inputs
to transcript assembly: neither can place a read across an intron, so StringTie /
Cufflinks built on their BAMs yield fragmented, mostly mono-exonic transcripts and
the lncRNA set silently degrades. `config/compatibility.yaml` therefore refuses the
combination unless the user overrides it, and steers them to HISAT2/STAR.

Legitimate uses here: aligning to a transcriptome/contig set, contamination
screening, or DNA-level work - not the RNA-seq assembly path.
"""
from __future__ import annotations
from pathlib import Path
from typing import List, Optional

from lnckit.utils.shell import RunResult
from lnckit.wrappers.base import BaseWrapper
from lnckit.wrappers.registry import register


@register("index", "bowtie2")
class Bowtie2BuildWrapper(BaseWrapper):
    requires = ("bowtie2-build",)
    tool_name = "bowtie2-build"

    def build_command(self, *, genome: Path, index_prefix: Path, **_) -> List[str]:
        return ["bowtie2-build", "--threads", self.threads,
                str(genome), str(index_prefix)]

    def expected_outputs(self, *, index_prefix: Path, **_) -> List[Path]:
        return [Path(f"{index_prefix}.1.bt2")]


@register("aligner", "bowtie2")
class Bowtie2AlignWrapper(BaseWrapper):
    requires = ("bowtie2",)
    tool_name = "bowtie2"

    def build_command(self, *, index_prefix: Path, sam_out: Optional[Path] = None,
                      fq1: Path, fq2: Optional[Path] = None, **_) -> List[str]:
        cmd = ["bowtie2", "-p", self.threads, "-x", str(index_prefix)]
        if fq2:
            cmd += ["-1", str(fq1), "-2", str(fq2)]
        else:
            cmd += ["-U", str(fq1)]
        # sam_out=None -> SAM to stdout (streamed into samtools sort); a path -> SAM file.
        if sam_out is not None:
            cmd += ["-S", str(sam_out)]
        return cmd

    def expected_outputs(self, *, sam_out: Optional[Path] = None, **_) -> List[Path]:
        return [Path(sam_out)] if sam_out is not None else []


@register("index", "bwa")
class BwaIndexWrapper(BaseWrapper):
    requires = ("bwa",)
    tool_name = "bwa-index"

    def build_command(self, *, genome: Path, index_prefix: Path, **_) -> List[str]:
        return ["bwa", "index", "-p", str(index_prefix), str(genome)]

    def expected_outputs(self, *, index_prefix: Path, **_) -> List[Path]:
        return [Path(f"{index_prefix}.bwt")]


@register("aligner", "bwa")
class BwaMemWrapper(BaseWrapper):
    """`bwa mem` writes SAM to stdout, so the SAM path is a stdout redirect."""

    requires = ("bwa",)
    tool_name = "bwa-mem"

    def build_command(self, *, index_prefix: Path, fq1: Path,
                      fq2: Optional[Path] = None, **_) -> List[str]:
        cmd = ["bwa", "mem", "-t", self.threads, str(index_prefix), str(fq1)]
        if fq2:
            cmd.append(str(fq2))
        return cmd

    def run(self, *, sam_out: Path, **kwargs) -> RunResult:
        return super().run(stdout_path=Path(sam_out), sam_out=sam_out, **kwargs)

    def expected_outputs(self, *, sam_out: Path, **_) -> List[Path]:
        return [Path(sam_out)]
