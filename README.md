<div align="center">

![Build Status](https://github.com/aghaeifar/SpinWalk/workflows/CMake/badge.svg)
![Lates Release](https://img.shields.io/github/v/release/aghaeifar/SpinWalk)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/aghaeifar/SpinWalk)
![GitHub top language](https://img.shields.io/github/languages/top/aghaeifar/SpinWalk)
![License](https://img.shields.io/github/license/aghaeifar/SpinWalk)

</div>

<p align="center">
<a href="https://github.com/aghaeifar/SpinWalk">
<img src="doc/img/logo.png" alt="MC/DC logo" width="160" height="180">
</a>
</p>

<p align="center">
<strong>Spins Random Walk Simulator</strong>
<br />
</p>

This program is designed to simulate the behavior of spins under a specific MR sequence within a microvascular network or digitally generated phantom. Breifly, the susceptibility variation between blood and tissue leads to local field inhomogeneity which accordingly can be used to generate an MR contrast. The program tries to perform a Monte-Carlo simulation over range of spins which are randomly or regularly distributed and diffuse in presence of user defined magnetic field. SpinWalk is a versatile Monte Carlo simulator. While its primary goal is to simulate BOLD fMRI, it can also be utilized for a wide range of diffusion applications.

Here are example plots obtained from SpinWalk where show BOLD sensitivity as a function of vessel size for Gradient Echo (GRE) and Spin Echo (SE) seqeuences.

![](./doc/img/gre_se.jpg)

Some [literature](#Literature) are provided as reference to get a better feeling of what are intended to get from this kind of simulations.

Simulator is written in C++ and utilizes CUDA to run in GPU. Therefore, it is possible to run a simulation within a short period of time provided that a good GPU presents. Simulator can detect all GPU cards, if there is more than one, and can distribute the task to run in parallel in multiple GPUs. This is helpful if you wan to run the simulator in HPC cluster with mutliple GPUs available in a node.


## How to build
You can either use the supplied Dockerfile for a consistent build that includes all dependencies, or install the dependencies separately and compile using CMake.

### Docker 
We suggest utilizing the provided Dockerfile, which automates the installation of all dependencies, as well as the cloning and building of the program. Download the [Dockerfile](./Dockerfile) to your current directory and then execute the following commands:

```
docker build -t spinwalk .
docker run --gpus all --rm -it --runtime=nvidia spinwalk bash
```
Execute the **spinwalk** command, and you'll encounter the help menu along with the list of available GPU(s) in the output.
### CMake

#### Dependencies

- A C++ compiler supprting C++ 17
- CUDA driver (*nvidia-smi* and *nvcc --version* must run in terminal)
- Boost libraries ([+](https://www.boost.org/))
- HDF5 Library ([+](https://www.hdfgroup.org/downloads/hdf5))
  
If you prefer to install the program without using Docker, follow these steps (tested in Ubuntu 22.04):

```
sudo apt-get update && apt-get install -y libboost-all-dev libhdf5-dev
git clone https://github.com/aghaeifar/SpinWalk.git
cd SpinWalk
cmake -B ./build
cmake --build ./build --config Release
```

## Quick test after build
After building the program, simply launch `spinwalk` in the terminal to view the help menu with all available options. At the end of the output, the detected GPUs and driver must be displayed.

## How to simulate

```
./spinwalk -c my_config.ini 
```

Several configuration files can be simulated sequentially:

```
./spinwalk -c config1.ini config2.ini ...
```


## Configuration files

Configruation file is a text based [ini file](https://en.wikipedia.org/wiki/INI_file) used to provide simulation parameters for simulator. Simulator can accept more than one configuration to simulation several configurations. A configuration file can inherit from another configuration file to avoid writing repetitive simulation parmeters. All the possible parameters are provided in [config_default.ini](./config/config_default.ini) with definition and expected unit.

## Input/Output file format
Spinwalk processes three distinct input files and produces a single output file. All files are stored on disk in _Hierarchical Data Format version 5_ ([HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format)) format. Reading and writing HDF5 file is natively supported by [MATLAB](https://www.mathworks.com/help/matlab/hdf5-files.html). Pythonic interface to work with HDF5 files is also avaialbe through [H5py](https://pypi.org/project/h5py/) package. Input files can be specified within the [FILES] section of the configuration file, while output filename is automatically generated and is stored in the output folder defined under the same [FILES] section in the configuration file.

### Inputs
#### Phantom

The simulation requires at least one phantom file. The phantom file comprises both an off-resonance map and a mask. The mask defines the geometry and types of various substrates. The inclusion of the off-resonance map is optional. When the off-resonance map is not included, it is assumed to be zero. The off-resonance must be normalized to the static magnetic field where is intended to be used for simulation (i.e., fieldmap must be calculated for 1T). It will be internally scaled to the given B0 paramater in configuration file. The phantom file is stored in HDF5 format and includes following datasets:

- `fov` : 3x1 single-precision floating-point array containing fov of 3D sample. Unit is meter. 
- `mask` : 3D 8-bit unsigned integer array indicating substrate types. The values within the mask are utilized to distinguish between various substrates (i.e., tissue or object types). Each substrate type must have its own diffusivity, T1 and T2 values defined in the configuration files.
- `fieldmap` : (optional) 3D single-precision floating-point array indicating off-resonance map. Unit is Tesla.
  
Size of mask and fieldmap must match.

Example Python script to write a phantom file can be like:
```python
import numpy as np
import h5py

N = 256;
fieldmap = np.random.rand(N, N, N).astype(np.float32)
mask = np.random.randint(0, 16, size=(N, N, N), dtype=np.uint8)
fov = np.array([1e-3, 1e-3, 1e-3], dtype=np.float32)

filename = './phantom.h5'
with h5py.File(filename, 'w') as h5f:
    h5f.create_dataset('mask', data=mask)
    h5f.create_dataset('fieldmap', data=fieldmap)
    h5f.create_dataset('fov', data=fov)
```

Example MATLAB script to write a phantom file can be like:
```matlab
N = 256;
fieldmap = single(rand(N, N, N)); 
mask = int8(randi([0, 15], N, N, N)); % MATLAB does not save uint8 in H5
fov = single([1e-3, 1e-3, 1e-3]); 

filename = './phantom.h5';

% Create HDF5 file and write datasets
h5create(filename, '/mask', size(mask), 'Datatype', 'int8');
h5write(filename, '/mask', permute(mask, 3:-1:1));

h5create(filename, '/fieldmap', size(fieldmap), 'Datatype', 'single');
h5write(filename, '/fieldmap', permute(fieldmap, 3:-1:1));

h5create(filename, '/fov', size(fov), 'Datatype', 'single');
h5write(filename, '/fov', fov);

```
**Note:** The memory layout used to store data in HDF5 is [row-major](https://en.wikipedia.org/wiki/Row-_and_column-major_order). As MATLAB follows a column-major standard, permuting dimensions is necessary when using MATLAB.

Multiple phantom files can be defined in the configuration file. Spinwalk will simulate them sequentially.

SpinWalk can also generate a digital phantom populated with cylinders or spheres. Details about the parameters are available in the help menu (execute *spinwalk -h*).

#### M0 and XYZ0 [optional]

M0 and XYZ0 are two additional inputs in configuration file which define starting magnization and initial spatial position of spins, respectively. These two are optional inputs, if not set or empty, spins will be positioned randomly with M0 = [0, 0, 1].

M0 and XYZ0 are of size *number of spins * 3* and stored as single-precision floating-point array in HDF5 file. Example MATLAB script to write M0 and XYZ0 file can be like:

```matlab
filename = 'XYZ0.h5';
h5create(filename, "/XYZ", size(XYZ0), 'Datatype','single');
h5write(filename, "/XYZ", permute(XYZ0, [2,1]));
filename = 'M0.h5';
h5create(filename, "/M", size(M0), 'Datatype','single');
h5write(filename, "/M", permute(M0, [2,1]));
```

unit for spatial position  in XYZ0 is meter.

### Outputs

The simulator generates a single output for each phantom. The output contains: 
- `M` : 4D single-precision floating-point array indicating the magnetization of spins at echo time(s). The dimentions are [scales, spins, echos, xyz (3)]
- `T` : 4D 8-bit unsigned integer array indicating the tissue or object type where spins are located at the time of the echo. The dimentions are [scales, spins, echos, 1]
- `XYZ` : 4D single-precision floating-point array indicating thethe spatial positions for either the entire random walk or just the final position. The dimentions are [scales, spins, steps, xyz (3)]
- `scales` : 1D double-precision floating-point array indicating scales used to scale the length/FoV of sample.
  
The folder to save outputs can be specified in the configuration file. the filename pattern is *{seqname}_{phantom_file_name}.h5*, where *seqname* and *phantom_file_name* are both defined by user in the configuration file.

Example MATLAB script to read output file:

```matlab
filename = 'output.h5';
M1 = h5read(filename, '/M');
M1 = permute(M1, 4:-1:1);
scales = h5read(filename, '/scales');
T = h5read(filename, '/T');  
T = permute(T, 4:-1:1);
```
## Literature

There are many nice papers published about simulation of BOLD signal in vessels network. A few are listed here for reference:

- Bieri O, Scheffler K. Effect of diffusion in inhomogeneous magnetic fields on balanced steady-state free precession. NMR Biomed. 2007 Feb;20(1):1-10. doi: 10.1002/nbm.1079. PMID: 16947639.
- Báez-Yánez MG, Ehses P, Mirkes C, Tsai PS, Kleinfeld D, Scheffler K. The impact of vessel size, orientation and intravascular contribution on the neurovascular fingerprint of BOLD bSSFP fMRI. Neuroimage. 2017 Dec;163:13-23. doi: 10.1016/j.neuroimage.2017.09.015. Epub 2017 Sep 8. PMID: 28890417; PMCID: PMC5857886.
- Boxerman JL, Hamberg LM, Rosen BR, Weisskoff RM. MR contrast due to intravascular magnetic susceptibility perturbations. Magn Reson Med. 1995 Oct;34(4):555-66. doi: 10.1002/mrm.1910340412. PMID: 8524024.
- Khajehim M, Nasiraei Moghaddam A. Investigating the spatial specificity of S2-SSFP fMRI: A Monte Carlo simulation approach. Magn Reson Imaging. 2017 Apr;37:282-289. doi: 10.1016/j.mri.2016.11.016. Epub 2016 Nov 24. PMID: 27890778.
- Weisskoff RM, Zuo CS, Boxerman JL, Rosen BR. Microscopic susceptibility variation and transverse relaxation: theory and experiment. Magn Reson Med. 1994 Jun;31(6):601-10. doi: 10.1002/mrm.1910310605. PMID: 8057812.
- Scheffler K, Engelmann J, Heule R. BOLD sensitivity and vessel size specificity along CPMG and GRASE echo trains. Magn Reson Med. 2021 Oct;86(4):2076-2083. doi: 10.1002/mrm.28871. Epub 2021 May 31. PMID: 34056746.

## Contributing

Contributions from the community to enhance and improve this project are welcome. For major changes, please open an issue first to discuss what you would like to change.