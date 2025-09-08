# Use a base image with a recent version of Ubuntu or Debian.
# Apptainer is well-supported on these distributions.
FROM ubuntu:22.04

# Install Apptainer and its dependencies.
# The `lsb-release` package is needed for the Apptainer installation script.
RUN apt-get update && \
    apt-get install -y \
    software-properties-common \
    lsb-release && \
    add-apt-repository -y ppa:apptainer/ppa && \
    apt-get update && \
    apt-get install -y apptainer && \
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