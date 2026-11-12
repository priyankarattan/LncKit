"""Infernal wrappers - Rfam covariance-model search for housekeeping ncRNAs.

Why this stage exists
---------------------
Transcript assembly does not only recover lncRNAs: the non-coding fraction is full
of rRNA, tRNA, snRNA, snoRNA and miRNA precursors. They survive every
coding-potential test *precisely because they genuinely are non-coding*, so
CPAT/CPC2 and protein-homology filtering cannot remove them. Only a structure-aware
search against Rfam separates them from real lncRNAs. Plant-LncPipe omits this
step, which is why its final set can still contain housekeeping ncRNA.

Two backends serve the `housekeeping` stage: Infernal `cmscan` (structure +
sequence - the sensitive option, used by default) and BLASTN against an ncRNA FASTA
(sequence only; see wrappers/blast.py).

`Rfam.cm` must be indexed by `cmpress` before `cmscan` can read it; the housekeeping
filter runs that automatically when the index is missing.
"""
from __future__ import annotations
from pathlib import Path
from typing import List, Optional

from lnckit.wrappers.base import BaseWrapper
from lnckit.wrappers.registry import register


@register("cmpress", "infernal")
class CmpressWrapper(BaseWrapper):
    """Index an Rfam.cm file (one-off; writes .i1f/.i1i/.i1m/.i1p siblings)."""

    requires = ("cmpress",)
    tool_name = "cmpress"

    def build_command(self, *, rfam_cm: Path, **_) -> List[str]:
        return ["cmpress", str(rfam_cm)]

    def expected_outputs(self, *, rfam_cm: Path, **_) -> List[Path]:
        return [Path(f"{rfam_cm}.i1m")]


@register("housekeeping", "infernal")
class CmscanWrapper(BaseWrapper):
    requires = ("cmscan",)
    tool_name = "cmscan"
    allow_empty_outputs = True   # a tblout carrying only comment lines is a valid result

    def build_command(self, *, fasta: Path, rfam_cm: Path, tblout: Path,
                      evalue: float = 1e-5, use_ga: bool = True,
                      clanin: Optional[Path] = None, **_) -> List[str]:
        cmd = ["cmscan", "--cpu", self.threads, "--rfam", "--nohmmonly",
               "--fmt", "2", "--tblout", str(tblout)]
        if use_ga:
            # Rfam's curated per-family gathering threshold - the cutoff Rfam itself
            # uses. cmscan treats --cut_ga and -E as mutually exclusive reporting
            # thresholds and exits with a usage error if both are passed, so this is
            # an either/or rather than a pair.
            cmd.append("--cut_ga")
        else:
            cmd += ["-E", evalue]
        if clanin:
            # Clan information lets cmscan mark redundant hits from the same Rfam
            # clan, so one locus is not counted as several families.
            cmd += ["--oclan", "--clanin", str(clanin)]
        cmd += [str(rfam_cm), str(fasta)]
        return cmd

    def expected_outputs(self, *, tblout: Path, **_) -> List[Path]:
        return [Path(tblout)]
