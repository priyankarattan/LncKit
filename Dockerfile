# LncKit container - every pipeline tool baked in. Build once, run anywhere with NO
# conda, NO compiling, NO per-machine install. Databases (UniProt/Rfam) and your data
# are bind-mounted at run time; they are NOT baked into the image.
#
# Build:  docker build -t lnckit:latest .
# Run:    docker run --rm -v /my/project:/work -v /my/db:/db lnckit:latest /work
FROM docker.io/condaforge/miniforge3:latest

SHELL ["/bin/bash", "-lc"]

# system build tools (CNCI/PLEK are compiled from source) + basics + a timezone
# (R warns loudly without one during the LncFinder install)
ENV TZ=UTC
RUN apt-get update && apt-get install -y --no-install-recommends \
        git curl unzip make g++ procps ca-certificates tzdata && \
    ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone && \
    rm -rf /var/lib/apt/lists/*

# Managed source tools (CNCI/PLEK) go to a FIXED baked path so they resolve for ANY
# user at run time (not $HOME-relative). Databases are mounted at /db.
ENV LNCKIT_HOME=/opt/lnckit_share
ENV PYTHONPATH=/opt/lnckit
ENV LNCKIT_DB_DIR=/db

WORKDIR /opt/lnckit
COPY . /opt/lnckit

# 1) build the 3 conda envs + FEELnc libdb fix + LncFinder (reuses setup.sh)
RUN bash setup.sh

# 2) build CNCI + PLEK from source so nothing installs at run time. Call the env's
# python directly (conda activate is unreliable in a non-interactive build) and put
# conda on PATH so the installer can create the Python-2 env CNCI/PLEK need.
RUN export PATH=/opt/conda/bin:$PATH && \
    /opt/conda/envs/lnckit/bin/python -c "from lnckit.core.installer import execute_plan, plan_install; [print(x) for x in execute_plan(plan_install(['cnci','plek']), yes=True)]" && \
    /opt/conda/bin/conda clean -afy

# run the pipeline in the directory passed as the first argument (your bind-mounted data)
ENTRYPOINT ["bash", "/opt/lnckit/lnckit.sh"]
