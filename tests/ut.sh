#!/bin/bash

# Get the directory where the script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Change to the script directory
pushd "$SCRIPT_DIR" > /dev/null

spinwalk -c ./../config/gre.ini  ./../config/se.ini ./../config/ssfp.ini ./../config/grase.ini ./../config/stimulated_echo.ini ./../config/gre_save_trajectory.ini 

