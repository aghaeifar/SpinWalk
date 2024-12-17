#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e
# Get a list of Docker images starting with "spinwalk"
images=$(docker images --format "{{.Repository}}:{{.Tag}}" | grep "^spinwalk")
TAR_FILE="/DATA2/spinwalk.tar"

for image in $images; do
    IMAGE_NAME=$image    
    APPTAINER_IMAGE="/DATA2/$image.sif"
    # Save the Docker image to a tar file
    echo "Saving Docker image to $TAR_FILE..."
    docker save -o $TAR_FILE $IMAGE_NAME
    # Load the Docker image tar file into Apptainer and convert it to a SIF image
    echo "Converting Docker tar file to Apptainer image..."
    apptainer build -F $APPTAINER_IMAGE docker-archive:$TAR_FILE
    echo "Apptainer image saved as $APPTAINER_IMAGE successfully."
done

