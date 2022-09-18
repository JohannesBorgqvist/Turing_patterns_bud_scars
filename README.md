



# Turing patterns on a sphere with holes, cell growth and local parameter changes
*Date:* 2021-09-01<br>
*Written by:* Johannes Borgqvist<br>


Description of the project

https://user-images.githubusercontent.com/77111216/190901847-791e2401-197a-4cc3-afd8-a9064b8fbf3a.mp4














## Working with branches
We have one main branch and two additional individual branches. All the work that should be done individually should be managed on the two individual branches called "Johannes\_branch" and "Philips\_branch" respectively. All the work that has been done in parallel should then be merged to the main branch. A good rule of thumb is that *no work is done on the main branch*, and we only merge to this branch. However, there are two important exceptions to this rule:

1. The README.md file should only be changed on the main branch,
2. The TASK_LIST.md file should only be changed on the main branch as well. 

Here is a good guide on *how to work with branches* ([see the following link](https://thenewstack.io/dont-mess-with-the-master-working-with-branches-in-git-and-github/)).



## A note on reproducibility
For the purpose of reproducibility, the project is entirely written in Python and the installation of all required packages has been managed through [anaconda](https://docs.anaconda.com/anaconda/install/index.html). Before presenting the instructions for installing all relevant packages using anaconda some information about the machine used to generate the results is presented. 

### Information about the laptop used to generate the results and write the code
The code has been developed and tested on a computer laptop with the following cpu information:

* Architecture:                    x86_64,
* CPU op-mode(s):                  32-bit, 64-bit,
* Byte Order:                      Little Endian,                                                                                                      
* Address sizes:                   39 bits physical, 48 bits virtual,                                                                                  
* CPU(s):                          8,
* On-line CPU(s) list:             0-7,
* Thread(s) per core:              2,
* Core(s) per socket:              4,
* Socket(s):                       1,
* NUMA node(s):                    1,
* Vendor ID:                       GenuineIntel,
* CPU family:                      6,
* Model:                           142,
* Model name:                      Intel(R) Core(TM) i7-10510U CPU @ 1.80GHz,
* Stepping:                        12,
* CPU MHz:                         800.057,
* CPU max MHz:                     4900,0000,
* CPU min MHz:                     400,0000,
* BogoMIPS:                        4599.93,
* Virtualisation:                  VT-x,
* L1d cache:                       128 KiB,
* L1i cache:                       128 KiB,
* L2 cache:                        1 MiB,
* L3 cache:                        8 MiB.

Some information about the *OS* of the machine is the following:

* Operating System: Ubuntu 21.04, 
* Kernel: Linux 5.11.0-31-generic,
* Architecture: x86-64.



### Installation of the packages using anaconda
In order to make the project as reproducible as possible, the coding has been done entirely in Python and the installation of all necessary packages is made possible by *anaconda*. The scripts associated with this project is enabled by three major platforms:
	
1. The [*FEniCS Project*](https://fenicsproject.org/) which is a "*popular open-source (LGPLv3) computing platform for solving partial differential equations (PDEs)*",
2. [*Gmsh*](https://gmsh.info/) which is a "*a three-dimensional finite element mesh generator with built-in pre- and post-processing facilities*",
3. [*ParaView*](https://www.paraview.org/) which is an "*open-source, multi-platform data analysis and visualization application*".

Both of the first two of these has a Python interface. Thus the scripts for generating the meshes using Gmsh and the scripts for conducting the FEM simulations using FEniCS are written in Python. For visualising the results, the graphical interface of ParaView has been used. 


Provided that anaconda has been succesfully installed, the easiest way to install all relevant packages is to install two different conda environments. These can be accessed through the provided yml-files called *gmsh\_latest\_version.yml* and *fenicsproject.yml*. These two environments are installed using the command<br>

*conda env create -f gmsh\_latest\_version.yml*<br>

and<br>
*conda env create -f fenicsproject.yml*<br>

respectively. The first environment must be activated in order to generate the meshes using the relevant Python scripts via Gmsh. The second environment must be activated in order to run the FEM based simulations using FEniCS. 

If you cannot install these to environments from the provided yml files, you can try to do it manually using the [*conda installation of FEniCS*](https://fenicsproject.org/download/). It is recommended to create and activate a conda environment in order to install and access all relevant parts of FEniCS. To create and activate a conda environment by the name of "*fenicsproject*", use the following two steps:

1. conda create -n fenicsproject -c conda-forge fenics
2. source activate fenicsproject

The first point installs all packages involved in the fenicsproject and stores them in a conda environment called "*fenicsproject*". The installation process will most likely take hours to complete. The second point activates the conda environment and when this environment is activated all important packages necessary to run FEM simulations using FEniCS are installed. 

In order to be able to generate the meshes and to visualise the results of the FEM-simulations, Gmsh and ParaView must be added to the conda environment "*fenicsproject*". This is done by typing the following four commands:

1. conda config --add channels conda-forge
2. conda config --set channel_priority strict
3. conda install --name fenicsproject gmsh python-gmsh
4. conda install --name fenicsproject paraview


When everythin works properly, it should be possible to run all scripts inside the conda environment "*fenicsproject*". This environment is activated using the command

* conda activate fenicsproject

and to exit the environment type

* conda deactivate

Lastly, as these platforms are under constant development the versions of the involved packages changes often. Therefore, the versions of the packages involved in this project have been documented in the Markdown-file "*VERSIONS\_OF\_PACKAGES.md*" in this repositry. 


Unfortunately, this approach does not install the latest version of Gmsh which is needed. So type the following command to create *another* environment in which the latest version of Gmsh is installed:

conda create -n gmsh\_latest\_version -c conda-forge python=3.9 gmsh=4.8 python-gmsh

Now, the environment "gmsh\_latest\_version" is entered when the meshes are created and "fenicsproject" is entered when we want to simulate the PDE system. 



### REALLY UGLY FIX FOR JOHANNES WORK COMPUTER

For some reason, "*conda activate*" does not work on this computer. So instead I need to run the following command: <br>

source /home/borgqvist/anaconda3/bin/activate gmsh\_latest\_version

It is also strange that "*conda deactivate*" seems to work, but apparently this command should be more or less equivalent

source /home/borgqvist/anaconda3/bin/deactivate gmsh\_latest\_version
