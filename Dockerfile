FROM ubuntu:18.04

####################################################################################
# Set environment variables
####################################################################################

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH
ENV JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/jre
ENV TZ=America/Denver

####################################################################################
# Install dependencies
####################################################################################

RUN apt-get update --fix-missing && \
  apt-get install -y wget curl git apt-transport-https software-properties-common && \
  apt-get update && \
  apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
  add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' && \
  apt-get update && \
  ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone && \
  apt-get -y --allow-unauthenticated install r-base r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev && \
  apt-get -y autoremove && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

####################################################################################
# Install R packages
####################################################################################

COPY install.R /
RUN Rscript /install.R

####################################################################################
# Install Python packages
####################################################################################

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh -O ~/miniconda.sh && \
  /bin/bash ~/miniconda.sh -b -p /opt/conda && \
  rm ~/miniconda.sh && \
  /opt/conda/bin/conda clean -tipsy && \
  ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
  echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
  echo "conda activate base" >> ~/.bashrc && \
  conda install numpy scikit-learn pandas tensorflow=1.11.0

####################################################################################
# Install Confounded within the image
####################################################################################

RUN cd /tmp && \
    git clone https://github.com/jdayton3/Confounded.git && \
    mkdir /confounded && \
    mv /tmp/Confounded/confounded /confounded && \
    mv /tmp/Confounded/data /confounded && \
    rm -rf /tmp/Confounded && \
    echo '#!/bin/bash\ncd /confounded\npython -m confounded "$@"' > /usr/bin/confounded && \
    chmod +x /usr/bin/confounded && \
    chmod 777 /confounded -R && \
    echo "Done importing Confounded code"

####################################################################################
# Copy files into the image
####################################################################################

COPY prepdata /prepdata
COPY adjust /adjust
COPY metrics /metrics
COPY figures /figures
COPY all.sh /

ENTRYPOINT /all.sh
