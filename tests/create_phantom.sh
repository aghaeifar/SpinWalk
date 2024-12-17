#!/bin/bash

# Get the directory where the script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Change to the script directory
pushd "$SCRIPT_DIR" > /dev/null

# Define the output directory
output_dir="./phantoms"
mkdir -p "$output_dir"

vol_frac=4
oxy_level_act=0.85
oxy_level_rest=0.78
dChi=0.00000011
orientation=90
resolution=600
fov=600
radius=8

# Define the output file
output_file_rest="${output_dir}/r${radius}_Y${oxy_level_rest}_vf${vol_frac}_ori${orientation}_fov${fov}_res${resolution}.h5"
# Call the command with the variable parameter and redirect the output
spinwalk phantom -c -r "$radius" -v "$vol_frac" -f "$fov" -z "$resolution" -d "$dChi" -y "$oxy_level_rest" -n "$orientation" -e 0 -o "$output_file_rest"

# Define the output file
output_file_act="${output_dir}/r${radius}_Y${oxy_level_act}_vf${vol_frac}_ori${orientation}_fov${fov}_res${resolution}.h5"
# Call the command with the variable parameter and redirect the output
spinwalk phantom -c -r "$radius" -v "$vol_frac" -f "$fov" -z "$resolution" -d "$dChi" -y "$oxy_level_act" -n "$orientation" -e 0 -o "$output_file_act"

# check phantoms are there
echo -e "\nList of existing phantoms in phantom folder..."
ls -l --block-size=M ./phantoms