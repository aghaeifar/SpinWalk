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


# Define the output file path
output_file="scales.ini"

# Write to the file
cat > "$output_file" <<EOL
[GENERAL]
PARENT_CONFIG = default_config.ini

[SIMULATION_PARAMETERS]
WHAT_TO_SCALE = 0
SCALE[0] = 0.0125
SCALE[1] = 0.0147
SCALE[2] = 0.0173
SCALE[3] = 0.0204
SCALE[4] = 0.0240
SCALE[5] = 0.0283
SCALE[6] = 0.0333
SCALE[7] = 0.0392
SCALE[8] = 0.0462
SCALE[9] = 0.0544
SCALE[10] = 0.0641
SCALE[11] = 0.0754
SCALE[12] = 0.0888
SCALE[13] = 0.1046
SCALE[14] = 0.1231
SCALE[15] = 0.1450
SCALE[16] = 0.1707
SCALE[17] = 0.2010
SCALE[18] = 0.2367
SCALE[19] = 0.2787
SCALE[20] = 0.3282
SCALE[21] = 0.3865
SCALE[22] = 0.4551
SCALE[23] = 0.5358
SCALE[24] = 0.6309
SCALE[25] = 0.7429
SCALE[26] = 0.8748
SCALE[27] = 1.0301
SCALE[28] = 1.2129
SCALE[29] = 1.4282
SCALE[30] = 1.6817
SCALE[31] = 1.9803
SCALE[32] = 2.3318
SCALE[33] = 2.7456
SCALE[34] = 3.2330
SCALE[35] = 3.8069
SCALE[36] = 4.4826
SCALE[37] = 5.2783
SCALE[38] = 6.2152
SCALE[39] = 7.3184
SCALE[40] = 8.6174
SCALE[41] = 10.1470
SCALE[42] = 11.9481
SCALE[43] = 14.0689
SCALE[44] = 16.5662
SCALE[45] = 19.5067
SCALE[46] = 22.9692
SCALE[47] = 27.0463
SCALE[48] = 31.8471
SCALE[49] = 37.5000
EOL
