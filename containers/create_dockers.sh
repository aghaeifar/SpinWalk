#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
pushd "$SCRIPT_DIR" > /dev/null

# docker system prune -a --volumes --force
docker system prune --force

CUDA_VERSIONS=("12.0.0" "12.6.3")

for VERSION in "${CUDA_VERSIONS[@]}"; do
    echo "Building Docker image for CUDA version $VERSION..."
    IMAGE_NAME="spinwalk_gpu:cuda${VERSION}"
    # Check if the image already exists
    if docker image inspect "$IMAGE_NAME" > /dev/null 2>&1; then
        echo "Image $IMAGE_NAME already exists. Deleting it..."
        docker rmi -f "$IMAGE_NAME"
    fi
    docker build  --no-cache -f ./Dockerfile -t "$IMAGE_NAME" --build-arg CUDA_VERSION=$VERSION ..
done


echo "Building CPU Docker image..."
IMAGE_NAME="spinwalk_cpu:ubuntu24.04"
if docker image inspect "$IMAGE_NAME" > /dev/null 2>&1; then
    echo "Image $IMAGE_NAME already exists. Deleting it..."
    docker rmi -f "$IMAGE_NAME"
fi
docker build --no-cache -f ./Dockerfile_CPU -t "$IMAGE_NAME" ..
