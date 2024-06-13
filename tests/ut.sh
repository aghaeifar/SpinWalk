#!/bin/bash

# Get the directory where the script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Change to the script directory
pushd "$SCRIPT_DIR" > /dev/null

spinwalk -c ./../config/gre.ini  ./../config/se.ini ./../config/ssfp.ini ./../config/grase.ini ./../config/stimulated_echo.ini ./../config/gre_save_trajectory.ini 

# spinwalk -c ./../config/dwi/dwi_baseline.ini ./../config/dwi/dwi_x_100.ini ./../config/dwi/dwi_x_200.ini ./../config/dwi/dwi_x_300.ini ./../config/dwi/dwi_x_400.ini ./../config/dwi/dwi_x_500.ini ./../config/dwi/dwi_x_600.ini ./../config/dwi/dwi_x_700.ini ./../config/dwi/dwi_x_800.ini ./../config/dwi/dwi_x_900.ini ./../config/dwi/dwi_x_1000.ini ./../config/dwi/dwi_x_1100.ini ./../config/dwi/dwi_x_1200.ini ./../config/dwi/dwi_x_1300.ini ./../config/dwi/dwi_x_1400.ini ./../config/dwi/dwi_x_1500.ini ./../config/dwi/dwi_x_1600.ini ./../config/dwi/dwi_x_1700.ini ./../config/dwi/dwi_x_1800.ini ./../config/dwi/dwi_x_1900.ini ./../config/dwi/dwi_x_2000.ini ./../config/dwi/dwi_x_2100.ini ./../config/dwi/dwi_x_2200.ini ./../config/dwi/dwi_x_2300.ini ./../config/dwi/dwi_x_2400.ini ./../config/dwi/dwi_x_2500.ini ./../config/dwi/dwi_x_2600.ini ./../config/dwi/dwi_x_2700.ini ./../config/dwi/dwi_x_2800.ini ./../config/dwi/dwi_x_2900.ini ./../config/dwi/dwi_x_3000.ini ./../config/dwi/dwi_x_3100.ini ./../config/dwi/dwi_x_3200.ini ./../config/dwi/dwi_x_3300.ini ./../config/dwi/dwi_x_3400.ini ./../config/dwi/dwi_x_3500.ini ./../config/dwi/dwi_x_3600.ini ./../config/dwi/dwi_x_3700.ini ./../config/dwi/dwi_x_3800.ini ./../config/dwi/dwi_x_3900.ini ./../config/dwi/dwi_x_4000.ini ./../config/dwi/dwi_x_4100.ini ./../config/dwi/dwi_x_4200.ini ./../config/dwi/dwi_x_4300.ini ./../config/dwi/dwi_x_4400.ini ./../config/dwi/dwi_x_4500.ini ./../config/dwi/dwi_x_4600.ini ./../config/dwi/dwi_x_4700.ini ./../config/dwi/dwi_x_4800.ini ./../config/dwi/dwi_x_4900.ini ./../config/dwi/dwi_x_5000.ini

