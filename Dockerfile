# FROM ubuntu:22.04
# FROM nvidia/cuda:12.0.1-runtime-ubuntu22.04
FROM nvidia/cuda:12.0.1-devel-ubuntu22.04

USER root

# needed for add-apt-repository
RUN apt-get update && apt install -y software-properties-common && apt-get update && apt-get clean && apt-get -y autoremove && rm -rf /var/lib/apt/lists/*  
# needed for gcc 11
RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test  && apt-get clean && apt-get -y autoremove && rm -rf /var/lib/apt/lists/*

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
    tmux \
	net-tools \
    iputils-ping \
	vim \
    git \
    libboost-all-dev \
    libhdf5-dev \
    python3-venv \
    && apt-get clean && apt-get -y autoremove && rm -rf /var/lib/apt/lists/* 

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

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir numpy torch scipy nibabel ipython jupyter matplotlib tqdm h5py pandas scikit-image scikit-learn seaborn virtualenv && rm -rf /root/.cache

# Clone the Git repository
COPY . /opt/SpinWalk/ 
# RUN git clone --depth 1 https://github.com/aghaeifar/SpinWalk.git /opt/SpinWalk
WORKDIR /opt/SpinWalk
RUN cmake -B ./build && cmake --build ./build --config Release
RUN cmake --install ./build


LABEL org.opencontainers.image.authors="Ali Aghaeifar"


# docker build -t spinwalk .
# docker run --gpus all --rm -it --runtime=nvidia spinwalk bash
