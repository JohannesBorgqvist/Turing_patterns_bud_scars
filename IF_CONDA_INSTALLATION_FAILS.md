*Date:* 2021-09-18<br>
*Written by:* Johannes Borgqvist<br>

If you cannot install the two relevant conda environments from the provided yml files, you can try to do it manually using the [*conda installation of FEniCS*](https://fenicsproject.org/download/). It is recommended to create and activate a conda environment in order to install and access all relevant parts of FEniCS. To create and activate a conda environment by the name of "*fenicsproject*", use the following two steps:

1. conda create -n fenicsproject -c conda-forge fenics
2. source activate fenicsproject

The first point installs all packages involved in the fenicsproject and stores them in a conda environment called "*fenicsproject*". The installation process will most likely take hours to complete. The second point activates the conda environment and when this environment is activated all important packages necessary to run FEM simulations using FEniCS are installed. 

In order to be able to generate the meshes and to visualise the results of the FEM-simulations, Gmsh and ParaView must be added to the conda environment "*fenicsproject*". This is done by typing the following four commands:

1. conda config --add channels conda-forge
2. conda config --set channel_priority strict
3. conda install --name fenicsproject gmsh python-gmsh
4. conda install --name fenicsproject paraview
. Lastly, as these platforms are under constant development the versions of the involved packages changes often. Therefore, the versions of the packages involved in this project have been documented in the Markdown-file "*VERSIONS\_OF\_PACKAGES.md*" in this repositry. 


Unfortunately, step 3 listed above does not install the latest version of Gmsh which is needed. So type the following command to create *another* environment in which the latest version of Gmsh is installed:

conda create -n gmsh\_latest\_version -c conda-forge python=3.9 gmsh=4.8 python-gmsh

This environment is also activate and deactivated exactly in the same manner as the environment "fenicsproject". The environment "gmsh\_latest\_version" is activated when the meshes are created and "fenicsproject" is activated when we want to simulate the PDE system. 
