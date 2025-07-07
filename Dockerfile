FROM bioconductor/bioconductor:RELEASE_3_19

####################################################################################
# Set environment variables
####################################################################################

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH
ENV JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/jre
ENV TZ=America/Denver

####################################################################################
# Env variables
####################################################################################

ENV TORCHINDUCTOR_CACHE_DIR=/tmp/torch_cache
# Using ARG to make the Conda path reusable and clean
ARG CONDA_DIR=/opt/conda
ENV PATH $CONDA_DIR/bin:$PATH
ENV CONDA_PREFIX=$CONDA_DIR

####################################################################################
# Install Miniforge3 FIRST to ensure Conda Python is primary
####################################################################################

RUN wget --quiet https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O ~/miniforge.sh && \
    /bin/bash ~/miniforge.sh -b -p $CONDA_DIR && \
    rm ~/miniforge.sh

# Explicitly tell R's reticulate package which Python to use.
ENV RETICULATE_PYTHON=$CONDA_DIR/bin/python

####################################################################################
# Create a user-writable directory for R packages and set ENV variable
####################################################################################

RUN mkdir -p /home/r_libs && \
    chmod -R 777 /home/r_libs
ENV R_LIBS_USER=/home/r_libs

####################################################################################
# Install R packages
####################################################################################

COPY install_main_packages.R /
RUN Rscript /install_main_packages.R
COPY install_annotation_packages.R /
RUN Rscript /install_annotation_packages.R
COPY install_adjuster_specific_packages.R /
RUN Rscript /install_adjuster_specific_packages.R

####################################################################################
# Install basic Python packages
####################################################################################

RUN conda install -y -c conda-forge numpy 'scikit-learn>=1.4' pandas 'tensorflow>=2.0' pytorch matplotlib

####################################################################################
# Install additional Python packages
####################################################################################

RUN conda install -y aif360 tabulate jaxtyping beartype aif360 pytables

####################################################################################
# Clone libraries
####################################################################################

RUN git clone https://github.com/edammer/TAMPOR.git /opt/TAMPOR
RUN git clone https://github.com/datapplab/AutoClass.git /opt/AutoClass

####################################################################################
# Pip
####################################################################################

RUN pip install 'aif360[FairAdapt]'

####################################################################################
# Install Confounded within the image
####################################################################################

#RUN cd /tmp && \
#    git clone https://github.com/jdayton3/Confounded.git && \
#    mkdir /confounded && \
#    mv /tmp/Confounded/confounded /confounded && \
#    mv /tmp/Confounded/data /confounded && \
#    rm -rf /tmp/Confounded && \
#    echo '#!/bin/bash\ncd /confounded\npython -m confounded "$@"' > /usr/bin/confounded && \
#    chmod +x /usr/bin/confounded && \
#    chmod 777 /confounded -R && \
#    echo "Done importing Confounded code"

