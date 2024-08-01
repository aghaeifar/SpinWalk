
FROM nvidia/cuda:12.5.1-devel-ubuntu22.04
USER root

# needed for add-apt-repository
RUN apt-get update && apt install -y software-properties-common && apt-get update && apt-get clean && apt-get -y autoremove && rm -rf /var/lib/apt/lists/*  

# Install the necessary dependencies for TBB and GCC 11
RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test  && \
    apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    libtbb-dev \
    gcc-11 \
    g++-11 \
	htop \
    tmux \
	net-tools \
    iputils-ping \
	vim \
    git \
    libboost-all-dev \
    libhdf5-dev \
    && apt-get clean && apt-get -y autoremove && rm -rf /var/lib/apt/lists/* 

# Set the default GCC version to 11
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 100 \
    && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-11 100

RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 11 && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-11 11

COPY . /opt/SpinWalk/ 
# RUN git clone --depth 1 https://github.com/aghaeifar/SpinWalk.git /opt/SpinWalk
WORKDIR /opt/SpinWalk
RUN cmake -B ./build && cmake --build ./build --config Release && cmake --install ./build

LABEL org.opencontainers.image.authors="Ali Aghaeifar"

# docker build -t spinwalk .
# docker run --gpus all --rm -it --runtime=nvidia spinwalk bash
