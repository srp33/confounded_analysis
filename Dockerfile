# Stage 1: Install system-level dependencies and ccache
FROM bioconductor/bioconductor:RELEASE_3_21 as builder-base

####################################################################################
# Set base environment variables
####################################################################################
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV TZ=America/Denver

####################################################################################
# Install system dependencies and ccache for compilation caching
####################################################################################
RUN apt-get update && \
    apt-get install -y --no-install-recommends wget gnupg software-properties-common git ccache && \
    wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null && \
    echo 'deb https://apt.kitware.com/ubuntu/ jammy main' | tee /etc/apt/sources.list.d/kitware.list >/dev/null && \
    apt-get update && \
    apt-get install -y --no-install-recommends cmake && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    cmake --version

####################################################################################
# Configure ccache to intercept and cache compiler calls
####################################################################################
ENV PATH /usr/lib/ccache:$PATH
ENV CC="ccache gcc"
ENV CXX="ccache g++"
# This directory will be mounted as a cache during the build
ENV CCACHE_DIR=/cache/ccache


# Stage 2: Build the Python environment in parallel
FROM builder-base as python-env

####################################################################################
# Set Conda environment variables
####################################################################################
ARG CONDA_DIR=/opt/conda
ENV PATH $CONDA_DIR/bin:$PATH
ENV CONDA_PREFIX=$CONDA_DIR

####################################################################################
# Install Miniforge3
####################################################################################
RUN wget --quiet https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O ~/miniforge.sh && \
    /bin/bash ~/miniforge.sh -b -p $CONDA_DIR && \
    rm ~/miniforge.sh

####################################################################################
# Install Python packages with persistent cache mounts for conda and pip
####################################################################################
RUN --mount=type=cache,id=conda-cache,target=/opt/conda/pkgs \
    --mount=type=cache,id=pip-cache,target=/root/.cache/pip \
    --mount=type=cache,id=ccache,target=/cache/ccache \
    conda install -y -c conda-forge numpy 'scikit-learn>=1.4' pandas 'tensorflow>=2.0' pytorch matplotlib && \
    conda install -y aif360 tabulate jaxtyping beartype pytables opentsne umap-learn geoparse seaborn && \
    conda install -y -c conda-forge psutil tqdm rich memory_profiler py-cpuinfo && \
    pip install gdown osfclient && \
    pip install 'aif360[FairAdapt]' && \
    conda clean --all -y

# Stage 3: Build the R environment in parallel
FROM builder-base as r-env

####################################################################################
# Set R library path
####################################################################################
ENV R_LIBS_USER=/home/r_libs
RUN mkdir -p $R_LIBS_USER && \
    chmod -R 777 $R_LIBS_USER

####################################################################################
# Install R packages using a single, consolidated script (Except Seurat)
####################################################################################
COPY install_packages.R /

####################################################################################
# Run the R installation with persistent cache mounts for pak and ccache
# FIX: Added --mount flags to cache R packages (pak) and compiled objects (ccache).
####################################################################################
RUN --mount=type=cache,id=r-pak-cache,target=/root/.cache \
    --mount=type=cache,id=ccache,target=/cache/ccache \
    Rscript /install_packages.R

# Final Stage: Assemble the final image
FROM builder-base

####################################################################################
# Copy artifacts from parallel stages
####################################################################################
COPY --from=python-env /opt/conda /opt/conda
COPY --from=r-env /home/r_libs /home/r_libs

####################################################################################
# Set all final environment variables
####################################################################################
ARG CONDA_DIR=/opt/conda
ENV PATH $CONDA_DIR/bin:$PATH
ENV CONDA_PREFIX=$CONDA_DIR
ENV RETICULATE_PYTHON=$CONDA_DIR/bin/python
ENV R_LIBS_USER=/home/r_libs
ENV TORCHINDUCTOR_CACHE_DIR=/tmp/torch_cache
ENV NUMBA_CACHE_DIR=/tmp/numba_cache

####################################################################################
# Create final directories and clone repositories
####################################################################################
RUN mkdir -p /tmp/numba_cache && chmod -R 777 /tmp/numba_cache
RUN git clone https://github.com/datapplab/AutoClass.git /opt/AutoClass
