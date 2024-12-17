#!/bin/bash

# Get the directory where the script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Change to the script directory
pushd "$SCRIPT_DIR" > /dev/null

spinwalk config -s GRE --TE 20000 -t 50 -p ./phantoms/r8_Y0.78_vf4_ori90_fov600_res600.h5 ./phantoms/r8_Y0.85_vf4_ori90_fov600_res600.h5 -o ./gre.ini
spinwalk config -s SE --TE 20000 -t 50 -p ./phantoms/r8_Y0.78_vf4_ori90_fov600_res600.h5 ./phantoms/r8_Y0.85_vf4_ori90_fov600_res600.h5 -o ./se.ini

# replace the default config file with the scales.ini as the parent config
sed -i 's/default_config.ini/scales.ini/g' gre.ini
sed -i 's/default_config.ini/scales.ini/g' se.ini
