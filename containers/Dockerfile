
# Stage 1: Build
ARG CUDA_VERSION=12.6.3
FROM nvidia/cuda:${CUDA_VERSION}-devel-ubuntu22.04 AS builder
USER root

# needed for add-apt-repository
RUN apt-get update && apt install -y software-properties-common --no-install-recommends && apt-get update && apt-get clean && apt-get -y autoremove && rm -rf /var/lib/apt/lists/*  

# Install the necessary dependencies for TBB and GCC 11
RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test  && \
    apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    libtbb-dev \
    gcc-11 \
    g++-11 \
    wget \ 
    libboost-all-dev \
    libhdf5-dev \
    && apt-get clean && apt-get -y autoremove && rm -rf /var/lib/apt/lists/* 

# Set the default GCC version to 11
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 100 \
    && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-11 100 \
    && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 11 \
    && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-11 11

# Set the CMake version
ARG CMAKE_VERSION=3.30.1
ARG CMAKE_DIR=cmake-${CMAKE_VERSION}-linux-x86_64

# Download and install CMake
RUN wget https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/${CMAKE_DIR}.tar.gz && \
    tar -zxvf ${CMAKE_DIR}.tar.gz && \
    mv ${CMAKE_DIR} /opt/cmake && \
    ln -s /opt/cmake/bin/* /usr/local/bin && \
    rm ${CMAKE_DIR}.tar.gz

COPY . /opt/SpinWalk/ 
WORKDIR /opt/SpinWalk
RUN rm -rf ./build && cmake -B ./build && cmake --build ./build --config Release && cmake --install ./build

# Stage 2: Runtime
FROM nvidia/cuda:${CUDA_VERSION}-runtime-ubuntu22.04
USER root
RUN apt-get update && apt install -y libboost-log1.74.0 libtbb-dev libhdf5-103 libgomp1 --no-install-recommends && apt-get clean && apt-get -y autoremove && rm -rf /var/lib/apt/lists/* 
COPY --from=builder /opt/SpinWalk/build/spinwalk /usr/local/bin/
LABEL org.opencontainers.image.authors="Ali Aghaeifar"

# docker build -t spinwalk .
# docker run --gpus all --rm -it --runtime=nvidia -v /DATA2:/mnt/DATA2 spinwalk_gpu:cuda12.6.3 bash
