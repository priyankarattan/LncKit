"""Auto-detect RNA-seq library strandedness from aligned reads.

Policy (safety by construction): auto-detect may only PROMOTE the library to a
stranded setting (RF or FR) when the evidence is overwhelming AND based on enough
informative reads. In every other case - weak or ambiguous signal, too few reads,
missing inputs, or ANY error - it returns None (unstranded). It never applies a
strand the data does not clearly support, so its worst case equals the safe
unstranded default, and an explicit --strandness from the user always wins (this
runs only when inputs.strandness == "auto").

Method (annotation-free): a splice-aware aligner tags each spliced read with the
transcript strand it came from - the XS:A: tag, inferred from the GT-AG splice
motif (HISAT2 emits it under --dta; STAR under --outSAMstrandField intronMotif).
Comparing each read's own mapped strand (the SAM FLAG) to that XS strand, using the
read1/read2 orientation rule from RSeQC infer_experiment, separates first-strand
(RF/dUTP) from second-strand (FR):

    read1 (or single-end): mapped-strand == XS -> FR vote ; != -> RF vote
    read2                : mapped-strand == XS -> RF vote ; != -> FR vote

An unstranded library splits the votes ~50/50; a stranded one is strongly one-sided.
Only spliced, primary, mapped reads that carry an XS tag are counted.
"""
from __future__ import annotations

import subprocess
from typing import Iterable, Optional, Tuple

from lnckit.core.logging_utils import get_logger

log = get_logger("lnckit.strandness")

#: need at least this many informative (spliced, XS-tagged) reads to make any call
MIN_INFORMATIVE = 20000
#: one orientation must reach this fraction of the votes to call RF or FR
MIN_FRACTION = 0.90
#: stop scanning after this many primary alignments (plenty for a confident call)
MAX_SCAN = 2_000_000


def _tally(sam_lines: Iterable[str]) -> Tuple[int, int, int]:
    """Count FR vs RF votes over spliced, XS-tagged primary reads.

    Pure function: takes SAM text lines, so the decision logic is unit-tested
    without needing an aligner or a BAM. Returns (n_informative, fr_votes, rf_votes).
    """
    n = fr = rf = 0
    for line in sam_lines:
        if not line or line[0] == "@":
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        try:
            flag = int(f[1])
        except ValueError:
            continue
        if flag & 0x4:                      # unmapped
            continue
        if "N" not in f[5]:                 # CIGAR has no splice -> no reliable XS
            continue
        xs = None
        for tag in f[11:]:
            if tag.startswith("XS:A:"):
                xs = tag[5:6]
                break
        if xs not in ("+", "-"):
            continue
        read_strand = "-" if (flag & 0x10) else "+"
        same = (read_strand == xs)
        if flag & 0x80:                     # read 2 of a pair: votes flip
            if same:
                rf += 1
            else:
                fr += 1
        else:                               # read 1, or single-end
            if same:
                fr += 1
            else:
                rf += 1
        n += 1
    return n, fr, rf


def decide(n_informative: int, fr_votes: int, rf_votes: int,
           min_informative: int = MIN_INFORMATIVE,
           min_fraction: float = MIN_FRACTION) -> Optional[str]:
    """Apply the confidence policy. Returns 'RF', 'FR', or None (unstranded)."""
    if n_informative < min_informative:
        return None
    p_fr = fr_votes / n_informative
    p_rf = rf_votes / n_informative
    if p_rf >= min_fraction:
        return "RF"
    if p_fr >= min_fraction:
        return "FR"
    return None


def _samtools_stream(bam, samtools: str = "samtools", max_scan: int = MAX_SCAN):
    """Yield SAM text lines for primary, mapped alignments (-F 0x104)."""
    proc = subprocess.Popen(
        [samtools, "view", "-F", "0x104", str(bam)],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    try:
        for i, line in enumerate(proc.stdout):
            if i >= max_scan:
                break
            yield line
    finally:
        try:
            proc.stdout.close()
        except Exception:
            pass
        try:
            proc.terminate()
        except Exception:
            pass


def detect(bam, samtools: str = "samtools") -> Tuple[Optional[str], str]:
    """Detect strandedness from a BAM. Returns (call, human_readable_reason).

    NEVER raises - any problem returns (None, reason) so the caller safely falls back
    to unstranded. `call` is 'RF', 'FR', or None.
    """
    try:
        n, fr, rf = _tally(_samtools_stream(bam, samtools))
    except Exception as e:                  # samtools missing, unreadable BAM, ...
        return None, f"detection error ({type(e).__name__}: {e}) - using unstranded"
    if n < MIN_INFORMATIVE:
        return None, (f"only {n} informative spliced reads (< {MIN_INFORMATIVE}) - "
                      f"using unstranded")
    p_fr, p_rf = fr / n, rf / n
    call = decide(n, fr, rf)
    if call == "RF":
        return "RF", (f"RF / first-strand (dUTP): {p_rf:.3f} of {n} informative reads "
                      f"are reverse-oriented")
    if call == "FR":
        return "FR", (f"FR / second-strand: {p_fr:.3f} of {n} informative reads are "
                      f"forward-oriented")
    return None, (f"ambiguous (FR {p_fr:.3f} / RF {p_rf:.3f} of {n} reads, neither >= "
                  f"{MIN_FRACTION}) - using unstranded")
