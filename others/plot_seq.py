#!/usr/bin/env python

import os
import sys
import configparser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def read_ini(file_path):
    config = configparser.ConfigParser()
    config.optionxform = str # Preserve case of keys
    try:
        config_parent_dict = {}
        config.read(file_path)
        if config.has_section('GENERAL'):
            if config.has_option('GENERAL', 'PARENT_CONFIG') and config.get('GENERAL', 'PARENT_CONFIG') != '':
                parent_config = config.get('GENERAL', 'PARENT_CONFIG')
                parent_config_path = os.path.join(os.path.dirname(file_path), parent_config)
                config_parent_dict = read_ini(parent_config_path)

        config_dict = {section: dict(config.items(section)) for section in config.sections()}
        config_combined_dict = {}
        all_keys = config_parent_dict.keys() | config_dict.keys()
        for key in all_keys:
            if key not in config_parent_dict:
                config_parent_dict[key] = {}
            if key not in config_dict:
                config_dict[key] = {}
            config_combined_dict[key] = {**config_parent_dict[key], **config_dict[key]}

    except Exception as e:
        print(f"Error reading INI file: {e}")
        sys.exit(1)

    return config_combined_dict

def main():
    if len(sys.argv) != 2:
        print("Usage: python plot_seq.py <file_path>")
        sys.exit(1)    
    file_path = sys.argv[1]

    if os.path.exists(file_path) is False :
        print(f"The file '{file_path}' does not exist.")
        sys.exit(1)
    
    config = read_ini(file_path)

    TR = np.float64(config['SCAN_PARAMETERS']['TR'])
    TE = np.fromstring(config['SCAN_PARAMETERS']["TE"], sep=' ')

    RF_T  = np.fromstring(config['SCAN_PARAMETERS']["RF_T"], sep=' ')
    RF_FA = np.fromstring(config['SCAN_PARAMETERS']["RF_FA"], sep=' ') 
    RF_PH = np.fromstring(config['SCAN_PARAMETERS']["RF_PH"], sep=' ') 
    
    GRADIENT_T = np.fromstring(config['SCAN_PARAMETERS']["GRADIENT_T"], sep=' ')
    GRADIENT_X = np.fromstring(config['SCAN_PARAMETERS']["GRADIENT_X"], sep=' ')
    GRADIENT_Y = np.fromstring(config['SCAN_PARAMETERS']["GRADIENT_Y"], sep=' ')
    GRADIENT_Z = np.fromstring(config['SCAN_PARAMETERS']["GRADIENT_Z"], sep=' ')

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6))
    for t, h in zip(RF_T, RF_FA):
        ax1.annotate('', xy=(t, h), xytext=(t, 0), arrowprops=dict(arrowstyle="->"))

    TE_h = np.max(RF_FA)*0.75
    for te in TE:
        ax1.add_patch(Rectangle((te - 50, 0), 100, TE_h, color='red'))

    ax1.set_xlim(-TR/10, TR)
    ax1.set_ylim(0, np.max(RF_FA)*1.1)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.set_ylabel('Flip Angle (deg)')

    ax2.plot(GRADIENT_T, GRADIENT_X, label='Gradient X')
    ax2.plot(GRADIENT_T, GRADIENT_Y, label='Gradient Y')
    ax2.plot(GRADIENT_T, GRADIENT_Z, label='Gradient Z')
    max_grad = np.max([np.max(GRADIENT_X), np.max(GRADIENT_Y), np.max(GRADIENT_Z)])
    min_grad = np.min([np.min(GRADIENT_X), np.min(GRADIENT_Y), np.min(GRADIENT_Z)])
    ax2.set_xlim(-TR/10, TR)
    ax2.set_ylim(min_grad*1.1, max_grad*1.1)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.legend()
    ax2.set_xlabel('Time (ms)')
    ax2.set_ylabel('Gradient Amplitude (mT/m)')

    plt.show()

if __name__ == "__main__":
    main()