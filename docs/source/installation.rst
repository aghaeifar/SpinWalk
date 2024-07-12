###############
Installation
###############


***************
How to build
***************

You can either use the supplied Dockerfile for a consistent build that includes all dependencies, or install the dependencies separately and compile using CMake.

Docker
===============

We suggest utilizing the provided Dockerfile, which automates the installation of all dependencies, as well as the cloning and building of the program. Download the Dockerfile to your current directory and then execute the following commands:


In almost any documentation you need to show examples like:


.. code-block:: bash

	docker build --no-cache -t spinwalk .
	docker run --gpus all --rm -it --runtime=nvidia spinwalk bash

CMake
===============

Dependencies
---------------
* A C++ compiler supprting C++ 20
* CUDA driver (nvidia-smi and nvcc --version must run in terminal)
* Boost libraries (+)
* HDF5 Library (+)
* Threading Building Blocks (TBB).

If you prefer to install the program without using Docker, follow these steps (tested in Ubuntu 22.04):

.. code-block:: bash

    sudo apt-get update && apt-get install -y libboost-all-dev libhdf5-dev libtbb-dev
    git clone https://github.com/aghaeifar/SpinWalk.git
    cd SpinWalk
    cmake -B ./build
    cmake --build ./build --config Release
    
**********************  
Quick test after build
**********************
After building the program, simply launch ``spinwalk`` in the terminal to view the help menu with all available options. At the end of the output, the detected GPUs and driver must be displayed.
