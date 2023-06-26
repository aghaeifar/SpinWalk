# Use a base image with MATLAB installed
# FROM mathworks/matlab:r2022b
# FROM ubuntu:22.04
# FROM nvidia/cuda:12.0.1-runtime-ubuntu22.04
 FROM nvidia/cuda:12.0.1-devel-ubuntu22.04

USER root

# needed for add-apt-repository
RUN apt-get update && apt install -y software-properties-common && apt-get update
# needed for gcc 11
RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test

# Install the necessary dependencies for TBB and GCC 11
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    libtbb-dev \
    gcc-11 \
    g++-11 \
    python3.10 \ 
    python3-pip \
    python-is-python3 \
	htop \
	net-tools \
	vim

# RUN apt-get install nvidia-cuda-toolkit

# Set the default GCC version to 11
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 100 \
    && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-11 100

RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 11 && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-11 11

# RUN  mkdir -p /mnt/host_home
# RUN  mkdir -p /mnt/host_tmp

RUN apt-get clean \
 && apt-get -y autoremove \
 && rm -rf /var/lib/apt/lists/*

LABEL org.opencontainers.image.authors="Ali Aghaeifar"



# Mount the host's home directory in the container


# Set the default command to start MATLAB
# CMD matlab
