"""Replicate-pooling stage - concatenate per-group TRIMMED FASTQs with `cat`.

Some designs pool biological replicates before assembly to maximise transcript-discovery
depth: e.g. all "5D_Bold-seeded" runs become ONE 5D_Bold-seeded sample, all
"10D_Bold-seeded" become ONE, etc. LncKit does this AFTER trimming - each replicate is
QC'd and trimmed individually, then the trimmed reads of every accession sharing a
`group` label in samples.txt are `cat`-ed into a single (pair of) file(s) per group under
`merged_reads/`.

One `cat` per group per mate: 4 groups of paired reads -> 4x2 = 8 `cat` calls, giving 4
merged samples. Gzipped FASTQ concatenation is valid (a concatenation of gzip streams is
itself a valid gzip stream). R1 files are cat'd in the SAME accession order as R2, so
mates stay in register.

Enabled with `inputs.merge_replicates: true` (or the wizard). No-op if no groups known.
"""
from __future__ import annotations

from typing import Dict, List

from lnckit.core.context import RunContext
from lnckit.core.logging_utils import get_logger
from lnckit.utils.shell import run as run_tool

log = get_logger("lnckit.merge")


def run_merge(ctx: RunContext, groups: Dict[str, List[str]],
              samples: List[Dict], paired: bool) -> List[Dict]:
    """Pool each group's TRIMMED samples into one merged sample.

    groups  : {group_label: [accession, ...]} from samples.txt
    samples : the trimmed per-run samples, each {"name": accession, "fq1":, "fq2":}
    Returns the merged sample list (one per group), each named after the group.
    """
    out = ctx.dirs["merged_reads"]
    resume = ctx.config.get("execution.resume", True)
    by_name = {s["name"]: s for s in samples}      # accession -> trimmed sample
    merged: List[Dict] = []

    for group, accs in groups.items():
        gname = str(group).strip().replace(" ", "_") or "group"
        members = [by_name[a] for a in accs if a in by_name]
        missing = [a for a in accs if a not in by_name]
        if missing:
            log.warning("group '%s': %d run(s) not found among trimmed samples: %s",
                        group, len(missing), ", ".join(missing))
        if not members:
            log.warning("group '%s': no trimmed runs found - skipping this group", group)
            continue

        if paired:
            r1 = [m["fq1"] for m in members if m.get("fq1")]
            r2 = [m["fq2"] for m in members if m.get("fq2")]
            if not r1 or len(r1) != len(r2):
                log.warning("group '%s': %d R1 / %d R2 file(s) - skipping merge",
                            group, len(r1), len(r2))
                continue
            m1, m2 = out / f"{gname}_1.fastq.gz", out / f"{gname}_2.fastq.gz"
            if not (m1.exists() and resume):
                run_tool(["cat", *map(str, r1)], tool="cat-merge", stdout_path=m1)
            if not (m2.exists() and resume):
                run_tool(["cat", *map(str, r2)], tool="cat-merge", stdout_path=m2)
            log.info("group '%s': merged %d run(s) -> %s + %s (from %s)",
                     group, len(r1), m1.name, m2.name,
                     ", ".join(m["name"] for m in members))
            merged.append({"name": gname, "fq1": m1, "fq2": m2})
        else:
            fs = [m["fq1"] for m in members if m.get("fq1")]
            if not fs:
                log.warning("group '%s': no reads - skipping", group)
                continue
            m = out / f"{gname}.fastq.gz"
            if not (m.exists() and resume):
                run_tool(["cat", *map(str, fs)], tool="cat-merge", stdout_path=m)
            log.info("group '%s': merged %d run(s) -> %s (from %s)",
                     group, len(fs), m.name, ", ".join(x["name"] for x in members))
            merged.append({"name": gname, "fq1": m, "fq2": None})

    if not merged:
        log.warning("merge_replicates was requested but no groups could be merged "
                    "(check the 'group' column in samples.txt)")
    return merged
