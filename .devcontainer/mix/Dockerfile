FROM quay.io/condaforge/miniforge3:23.3.1-1

# add conda-forge and bioconda channels to the top of the priority list and enable mamba solver
RUN /opt/conda/bin/conda config --add channels defaults && \
    /opt/conda/bin/conda config --add channels bioconda && \
    /opt/conda/bin/conda config --add channels conda-forge
    #  channel priority strict seems to mess up package resolution
    # /opt/conda/bin/conda config --set channel_priority strict

# Copy environment.yml (if found) to a temp location so we update the environment. Also
# copy "noop.txt" so the COPY instruction does not fail if no environment.yml exists.
COPY environment.yml* .devcontainer/noop.txt /tmp/conda-tmp/
RUN if [ -f "/tmp/conda-tmp/environment.yml" ]; then umask 0002 \
    && /opt/conda/bin/mamba env update -n base -f /tmp/conda-tmp/environment.yml; fi \
    && rm -rf /tmp/conda-tmp

RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends \
    nano \
    micro \
    vim \
    less \
    git \
    screen \
    tmux \
    wget \
    curl \
    colordiff \
    htop \
    tree \
    graphviz
