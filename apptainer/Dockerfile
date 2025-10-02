# Use a base image with a recent version of Ubuntu or Debian.
# Apptainer is well-supported on these distributions.
FROM ubuntu:24.04

# Install Apptainer and its dependencies using the official installation method
RUN apt-get update && \
    apt-get install -y \
    wget \
    ca-certificates && \
    wget -O apptainer.deb https://github.com/apptainer/apptainer/releases/download/v1.3.4/apptainer_1.3.4_amd64.deb && \
    apt-get install -y ./apptainer.deb && \
    rm apptainer.deb && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Copy the Apptainer definition file and the R installation script into the image.
# COPY apptainer.def /apptainer.def
# COPY install_packages.R /install_packages.R

# Set a working directory inside the container.
WORKDIR /

# Command to build the Apptainer image.
# This will execute when the container is run.
# The 'fakeroot' option is used to build a non-privileged SIF image.
CMD ["apptainer", "build", "--fakeroot", "my_app.sif", "apptainer.def"]