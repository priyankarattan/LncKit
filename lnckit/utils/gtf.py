"""Lightweight GTF helpers: transcript IDs, per-transcript exon counts, ID subsetting."""
from __future__ import annotations
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Set

_TID = re.compile(r'transcript_id "([^"]+)"')
_TPM = re.compile(r'TPM "([^"]+)"')
_FPKM = re.compile(r'FPKM "([^"]+)"')


def transcript_abundance(gtf: Path, metric: str = "tpm") -> Dict[str, float]:
    """transcript_id -> TPM or FPKM from a StringTie `-e` abundance GTF.

    StringTie writes both `TPM "..."` and `FPKM "..."` on each transcript line; `metric`
    selects which one (default TPM). FPKM = fragments per kilobase per million.
    """
    pat = _FPKM if str(metric).lower() == "fpkm" else _TPM
    out: Dict[str, float] = {}
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9 or parts[2] != "transcript":
                continue
            tid, val = _TID.search(line), pat.search(line)
            if tid and val:
                try:
                    out[tid.group(1)] = float(val.group(1))
                except ValueError:
                    continue
    return out


def transcript_tpm(gtf: Path) -> Dict[str, float]:
    """transcript_id -> TPM from a StringTie `-e` abundance GTF (back-compat)."""
    return transcript_abundance(gtf, "tpm")


def transcript_ids(gtf: Path) -> Set[str]:
    ids = set()
    with open(gtf) as fh:
        for line in fh:
            m = _TID.search(line)
            if m:
                ids.add(m.group(1))
    return ids


def exon_counts(gtf: Path) -> Dict[str, int]:
    counts: Dict[str, int] = defaultdict(int)
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) > 2 and parts[2] == "exon":
                m = _TID.search(line)
                if m:
                    counts[m.group(1)] += 1
    return dict(counts)


def subset_gtf_by_ids(gtf: Path, keep_ids: Set[str], out: Path) -> Path:
    with open(gtf) as fin, open(out, "w") as fout:
        for line in fin:
            m = _TID.search(line)
            if m and m.group(1) in keep_ids:
                fout.write(line)
    return out


def has_empty_transcript_ids(gtf: Path) -> bool:
    """True if any feature row carries an empty `transcript_id ""`.

    NCBI RefSeq GTFs write their `gene` rows this way; StringTie's -G reader aborts on
    them with 'Error: no valid ID found for GFF record'. Detecting this lets the engine
    substitute a cleaned copy only when it is actually needed.
    """
    with open(gtf) as fh:
        for line in fh:
            if not line.startswith("#") and 'transcript_id ""' in line:
                return True
    return False


#: The only strand values StringTie's -G reader accepts. NCBI marks trans-spliced
#: organellar genes (chloroplast rps12, some tRNAs) with '?' (undetermined strand),
#: and StringTie aborts on them with "Error parsing strand (?) from GFF line".
_VALID_STRANDS = {"+", "-", "."}


def _strand_of(line: str):
    """The strand (column 7) of a GTF feature row, or None for comments/short rows."""
    if line.startswith("#"):
        return None
    parts = line.split("\t")
    if len(parts) < 8:
        return None
    return parts[6]


def has_invalid_strand(gtf: Path) -> bool:
    """True if any feature row has a strand other than +, -, or . (e.g. '?').

    Complements has_empty_transcript_ids(): either condition means the annotation
    needs a cleaned copy before StringTie's -G reader will accept it.
    """
    with open(gtf) as fh:
        for line in fh:
            s = _strand_of(line)
            if s is not None and s not in _VALID_STRANDS:
                return True
    return False


def sanitize_gtf(gtf: Path, out: Path) -> Dict[str, int]:
    """Write `out` = `gtf` minus the rows StringTie's -G reader rejects.

    Two kinds of rows are removed:
      * empty/missing transcript_id - the NCBI `gene` records (StringTie: 'no valid ID
        found for GFF record'); they carry nothing StringTie needs.
      * every row of a transcript whose strand is undetermined ('?') or otherwise not
        +/-/. - trans-spliced organellar genes such as chloroplast rps12 (StringTie:
        'Error parsing strand (?)'). The WHOLE transcript is dropped, not just the
        offending line, so no orphan exons are left referencing a missing transcript.

    Returns a breakdown: {'no_tid': n, 'bad_strand': n, 'total': n}.
    """
    # Pass 1: collect transcript_ids that carry at least one invalid-strand row.
    bad_tids: Set[str] = set()
    with open(gtf) as fh:
        for line in fh:
            s = _strand_of(line)
            if s is not None and s not in _VALID_STRANDS:
                m = _TID.search(line)
                if m:
                    bad_tids.add(m.group(1))

    # Pass 2: write the survivors, tallying why each dropped row was removed.
    counts = {"no_tid": 0, "bad_strand": 0, "total": 0}
    with open(gtf) as fin, open(out, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            m = _TID.search(line)
            if not m:                                   # empty/missing transcript_id
                counts["no_tid"] += 1
            elif m.group(1) in bad_tids:                # trans-spliced / '?'-strand
                counts["bad_strand"] += 1
            else:
                fout.write(line)
                continue
            counts["total"] += 1
    return counts
