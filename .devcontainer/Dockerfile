FROM quay.io/condaforge/miniforge3:23.11.0-0

# Or your actual UID, GID on Linux if not the default 1000
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID="${USER_UID}"

# add conda-forge and bioconda channels to the top of the priority list and enable mamba solver
RUN /opt/conda/bin/conda config --add channels bioconda \
    && /opt/conda/bin/conda config --add channels conda-forge
#  channel priority strict seems to mess up package resolution
# /opt/conda/bin/conda config --set channel_priority strict

# Copy environment.yml (if found) to a temp location so we update the environment. Also
# copy "noop.txt" so the COPY instruction does not fail if no environment.yml exists.
# This updates the base environment instead of creating a new named virtual environment (mamba env create -f env.yml)
COPY environment.yml* .devcontainer/noop.txt /tmp/conda-tmp/
RUN if [ -f "/tmp/conda-tmp/environment.yml" ]; then umask 0002 \
    && /opt/conda/bin/mamba env update -n base -f /tmp/conda-tmp/environment.yml; fi \
    && rm -rf /tmp/conda-tmp

# Install additional packages.
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends \
    nano \
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

# Create a non-root user - see https://aka.ms/vscode-remote/containers/non-root-user,
# https://github.com/microsoft/vscode-dev-containers/issues/108,
# and https://github.com/mamba-org/micromamba-docker/blob/521003c092aca4be4800220a883ebebf5bacbd6d/Dockerfile#L44
RUN groupadd --gid 1000 "${USERNAME}" \
    && useradd -s /bin/bash --uid "${USER_UID}" --gid "${USER_GID}" -m "${USERNAME}" \
    # Add sudo support
    && apt-get install -y sudo \
    && echo "${USERNAME}" ALL=\(root\) NOPASSWD:ALL > "/etc/sudoers.d/${USERNAME}" \
    && chmod 0440 "/etc/sudoers.d/${USERNAME}" \
    # Transfer mamba ownership to user
    && chmod -R a+rwx /opt/conda \
    # Clean up
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*
