"""Genome-indexing stage (HISAT2 / STAR / Bowtie2 / BWA).

The backend is whatever `tools.aligner` names: each aligner registers its own index
wrapper under the `index` stage, so this module never branches on tool names except
where the output *shape* differs - STAR produces an index directory, the others a
file prefix.
"""
from __future__ import annotations

from lnckit.core.context import RunContext
from lnckit.core.exceptions import ConfigError
from lnckit.core.logging_utils import get_logger
from lnckit.wrappers import get_wrapper

log = get_logger("lnckit.index")


def run_index(ctx: RunContext):
    backend = ctx.config.tools.get("aligner", "hisat2")
    genome = ctx.config.get("reference.genome_fasta")
    if not genome:
        raise ConfigError("reference.genome_fasta is required to build the genome index.")

    idx_dir = ctx.dirs["index"]
    resume = ctx.config.get("execution.resume", True)
    wrapper = get_wrapper("index", backend)(threads=ctx.config.threads, log_dir=ctx.logs)

    if backend == "star":
        # STAR signals a finished index with SAindex; its presence means resume can skip.
        if resume and (idx_dir / "SAindex").exists():
            log.info("STAR index already present - skipping indexing")
        else:
            wrapper.run(genome=genome, index_dir=idx_dir,
                        annotation=ctx.config.get("reference.annotation_gtf"))
        return idx_dir

    prefix = idx_dir / "genome"
    # hisat2/bowtie2/bwa all write files that start with the prefix basename
    # (genome.*.ht2 / .bt2 / .bwt ...); any such file means the index is built.
    if resume and any(idx_dir.glob(f"{prefix.name}.*")):
        log.info("%s index already present - skipping indexing", backend)
    else:
        wrapper.run(genome=genome, index_prefix=prefix)
    return prefix
