FROM bioconductor/bioconductor:RELEASE_3_19

####################################################################################
# Set environment variables
####################################################################################

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH
ENV JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/jre
ENV TZ=America/Denver

####################################################################################
# Install R packages
####################################################################################

COPY install_*.R /
RUN Rscript /install_main_packages.R

####################################################################################
# Install Miniforge3 (which includes Mamba by default)
####################################################################################
RUN wget --quiet https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O ~/miniforge.sh && \
    /bin/bash ~/miniforge.sh -b -p /opt/conda && \
    rm ~/miniforge.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# RUN apt-get update --fix-missing && \
#  apt-get install -y wget curl git parallel apt-transport-https software-properties-common && \
#  apt-get update && \
#  apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
#  add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/' && \
#  apt-get update && \
#  ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone && \
#  apt-get -y --allow-unauthenticated install r-base r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev && \
#  apt-get -y autoremove && \
#  apt-get clean && \
#  rm -rf /var/lib/apt/lists/*

####################################################################################
# Install basic Python packages
####################################################################################

RUN conda install -y -c conda-forge numpy 'scikit-learn>=1.4' pandas 'tensorflow>=2.0' pytorch matplotlib

####################################################################################
# Install additional Python packages
####################################################################################

RUN conda install -y -c conda-forge tabulate
# RUN pip3 install numpy scikit-learn pandas 
# tensorflow=1.11.0

####################################################################################
# Install additional R packages
####################################################################################

RUN Rscript /install_annotation_packages.R
RUN Rscript /install_adjuster_specific_packages.R

####################################################################################
# Clone libraries
####################################################################################

RUN git clone https://github.com/edammer/TAMPOR.git /opt/TAMPOR
RUN git clone https://github.com/datapplab/AutoClass.git /opt/AutoClass

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

