"""Minimal, dependency-free FASTA utilities + a pure-Python ORF calculator.

The pure-Python longest-ORF function is a fallback so the ORF filter still works
in environments where no external ORF caller is installed (e.g. CI/tests). In
production runs the external getorf/TransDecoder/ORFfinder backends are used.
"""
from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterator, Tuple

_STOP = {"TAA", "TAG", "TGA"}
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def fasta_lengths(path: Path) -> Dict[str, int]:
    return {name: len(seq) for name, seq in read_fasta(path)}


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def subset_fasta(path: Path, keep: set, out: Path) -> Path:
    """Write only the named records to `out`.

    Each filter narrows the candidate set, so the next tool should only see the
    survivors: running CPAT/DIAMOND/cmscan over the full assembly instead of the
    current candidates is the difference between minutes and hours on a real
    transcriptome.
    """
    with open(out, "w") as fh:
        for name, seq in read_fasta(path):
            if name in keep:
                fh.write(f">{name}\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")
    return Path(out)


def longest_orf_length(seq: str, strand: str = "both") -> int:
    """Return length (in nt) of the longest ORF (ATG..STOP) across frames.

    Considers all three reading frames on the forward strand and, unless
    strand=='forward', the three frames on the reverse strand.
    """
    seq = seq.upper().replace("U", "T")
    best = 0
    strands = [seq] if strand == "forward" else [seq, revcomp(seq)]
    for s in strands:
        for frame in range(3):
            i = frame
            n = len(s)
            while i < n - 2:
                if s[i:i + 3] == "ATG":
                    j = i
                    while j < n - 2:
                        if s[j:j + 3] in _STOP:
                            best = max(best, j + 3 - i)
                            break
                        j += 3
                    i = j + 3
                else:
                    i += 3
    return best
