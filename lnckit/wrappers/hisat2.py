"""HISAT2 wrapper - splice-aware alignment (recommended for StringTie/Cufflinks)."""
from __future__ import annotations
from pathlib import Path
from typing import List, Optional
from lnckit.wrappers.base import BaseWrapper
from lnckit.wrappers.registry import register


@register("index", "hisat2")
class Hisat2BuildWrapper(BaseWrapper):
    requires = ("hisat2-build",)
    tool_name = "hisat2-build"

    def build_command(self, *, genome: Path, index_prefix: Path, **_) -> List[str]:
        return ["hisat2-build", "-p", self.threads, str(genome), str(index_prefix)]

    def expected_outputs(self, *, index_prefix: Path, **_) -> List[Path]:
        return [Path(f"{index_prefix}.1.ht2")]


@register("aligner", "hisat2")
class Hisat2AlignWrapper(BaseWrapper):
    requires = ("hisat2",)
    tool_name = "hisat2"

    def build_command(self, *, index_prefix: Path, sam_out: Optional[Path] = None,
                      fq1: Path, fq2: Optional[Path] = None,
                      strandness: Optional[str] = None, dta: bool = True,
                      **_) -> List[str]:
        # --dta ("downstream transcriptome assembly"): recommended by StringTie - without it
        # HISAT2 uses short anchors on novel splice sites and StringTie assembles more
        # fragmented transcripts. ON by default. Some published pipelines (e.g. Plant-LncPipe)
        # omit it; set tool_params.hisat2.dta=false to reproduce those runs exactly.
        cmd = ["hisat2", "--new-summary", "-p", self.threads, "-x", str(index_prefix)]
        if dta:
            cmd.insert(1, "--dta")
        if strandness:
            cmd += ["--rna-strandness", strandness]
        if fq2:
            cmd += ["-1", str(fq1), "-2", str(fq2)]
        else:
            cmd += ["-U", str(fq1)]
        # sam_out=None -> write SAM to stdout so it can be streamed straight into
        # `samtools sort` (no multi-GB SAM on disk). A path -> write that SAM file.
        if sam_out is not None:
            cmd += ["-S", str(sam_out)]
        return cmd

    def expected_outputs(self, *, sam_out: Optional[Path] = None, **_) -> List[Path]:
        return [Path(sam_out)] if sam_out is not None else []
