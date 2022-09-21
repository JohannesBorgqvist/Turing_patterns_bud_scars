# A hole in Turing’s theory: pattern formation on the sphere with budscars
*Date:* 2021-09-20<br>
*Written by:* Johannes Borgqvist<br>
Welcome to the github repositry associated with the article "*A hole in Turing’s theory: pattern formation on the sphere with budscars*" (**Add link to the future article**). The starting point of this project is the Schnakenberg model which has previously been analysed on the sphere [1]. If we set the initial conditions for both species of this model to small perturbations around the respective steady state concentrations, then patterns in the concentration profile will be formed according to the phenomena called *diffusion-driven instability* initially proposed by Alan M. Turing. In the case of the sphere, these patterns will correspond to the formation of *poles* which are defined as high-concentration regions of the active component (see the simulation below and the figure below).  



https://user-images.githubusercontent.com/77111216/190901847-791e2401-197a-4cc3-afd8-a9064b8fbf3a.mp4


![Pattern of the Schnakenberg model](./Figures/Schnakenberg_pattern_formation.png)


Now, in this project we are interested in adding a single hole in the spatial domain, i.e. the unit sphere, in order to see how this change in the domain affects the resulting patterns. To be able to generate all the relevant results of the article, you can run the bash-script called "*run\_all.sh*" (for instructions on how to do this, see the last section of this document). Alternatively, you can follow the "*STEP\_BY\_STEP\_GUIDE.md*" in the Code folder. Before you run the run\_all.sh bash-script or navigate to the code folder, you need to install all relevant packages, and how to achieve this is described in the section of this document that is entitled "*Installation of relevant packages and libraries using anaconda*".<br> 

Before this is described, I wish to say a word or two about my own philosophy on the notion of *reproducibility* which permeates this project.  




# A note on reproducibility
For the purpose of reproducibility, the project is entirely written in Python and the installation of all required packages is managed by [anaconda](https://docs.anaconda.com/anaconda/install/index.html). The version of Python that has been used is 3.8.3, and mainly two packages are frequently used:

1. *numpy*, version 1.18.5,
2. *matplotlib*, version 3.2.2.

In addition to these packages, there are numerous rather large platforms that are required, and these are installed using anaconda which will be described in the next section of this document. 

As I said before, this repositry is based on the idea of *reproducibility* which basically means that anyone should be able to generate the results presented in the article on their own. To promote this idea, several concrete steps have been taken:

1. The project is completely open-source and written in open-source languages (mainly Python).
2. The repositry has a clear folder structure where the folders are called Code, Figures, Meshes and Output.
3. The repositry is well-documented with numerous README.md files in each folder containing clear instructions on, for example, how to run the scripts as well as installing relevant packages.
4. The installation of the numerous platforms required for generating the results is simplified by the usage of conda (see the next section for more details).
5. The relevant functions carrying out specific tasks are divided into numerous scripts with clear names. For example, all scripts related to the *finite element method (FEM)* are stored in the script "*toolbox\_FEM\_simulations\_Schnakenberg\_sphere\_with\_holes.py*".
6. The code itself contains a lot of comments to clarify what is being done in a specific part of any particular script. 
7. We wrote a bash-script called "*run\_all.sh*" which runs all important scripts in succession and generate all results.


Before presenting the instructions for installing all relevant packages using anaconda, some information about the machine used to generate the results is presented. 

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



# Installation of relevant packages and libraries using anaconda

In order to make the project as reproducible as possible, the coding has been done entirely in Python and the installation of all necessary packages is made possible by the package [*conda*](https://anaconda.org/anaconda/conda) which is part of *anaconda*. The scripts associated with this project is enabled by four major platforms:
	
1. The [*FEniCS Project*](https://fenicsproject.org/) which is a "*popular open-source (LGPLv3) computing platform for solving partial differential equations (PDEs)*",
2. [*Gmsh*](https://gmsh.info/) which is a "*a three-dimensional finite element mesh generator with built-in pre- and post-processing facilities*",
3. [*scikit-learn*](https://scikit-learn.org/stable/) which is a library that enables "*machine Learning in Python*",
4. [*ParaView*](https://www.paraview.org/) which is an "*open-source, multi-platform data analysis and visualization application*".


The first of these two have Python interfaces, and the third is itself a library in Python. Thus, all scripts which involve the scripts for generating the meshes using Gmsh, the scripts for conducting the FEM simulations using FEniCS and the scripts for analysing the data using scikit-learn are written in Python. In order to visualise the results, the graphical interface ParaView has been used as well as [*PGFPlots*](http://pgfplots.sourceforge.net/) in LaTeX. In this project, we have not included the scripts for generating the LaTeX based plots, but instead we generated the same plots in Python using matplotlib. A diagram of how all the different parts of the project are connected is presented below.  

![Work flow](./Figures/diagram.png)


Provided that anaconda has been succesfully installed, the easiest way to install all relevant packages is to install two different conda environments. These can be accessed through the provided yml-files called *gmsh\_latest\_version.yml* and *fenicsproject.yml*. These two environments are installed using the command<br>

*conda env create -f gmsh\_latest\_version.yml*<br>

and<br>
*conda env create -f fenicsproject.yml*<br>

respectively. The first environment must be activated in order to generate the meshes using the relevant Python scripts that depend on Gmsh. The second environment must be activated in order to run the FEM based simulations using FEniCS. 



When everything works properly, it should be possible to run all scripts inside these two conda environments. Both environments are activated and deactivated in the same way. For example, the environment  "*fenicsproject*" is activated using the command

* *conda activate fenicsproject*

and to exit this environment type

* *conda deactivate*<br>

. On certain computers, you might have to replace the command "*conda activate <env\_name>*" with "*source activate <env\_name>*". As we said previously, the environment "*gmsh\_latest\_version*" is activated and deactivated in the same way.<br> 


If, God forbid, the installation of all relevant platforms based on these yml-files fails for some reason, all involved libraries must be installed manually. Some guidance on how this can be done is presented in the file "*IF\_CONDA\_INSTALLATION\_FAILS.md*".<br>

On the other hand, if you have managed to install the two conda environments "*fenicsproject*" and "*gmsh\_latest\_version*" you are now ready to start running simulations. 

# Running all of the scripts

The easiest way to run all relevant script is to use the bash-script called *run\_all.sh*. To run this script, begin by typing:<br>

"*chmod +x run\_all.sh*"<br>

which gives you permission to run all scripts in succession. You might need to add a "sudo" in front of chmod to make this work. After this, you need to change line 14 in the script "./run\_all.sh" from:<br>

*source /home/johannes/anaconda3/etc/profile.d/conda.sh*<br>

to the correct path where anaconda is installed on your computer. When you have entered the correct path, you can run all scripts in the following way:<br>

"*./run\_all.sh*"<br>

and this will generate all meshes, run the FEM-simulations as well as analyse the generated data. 



If you are interested in the individual steps, you can navigate to the Code-folder and follow the detailed instructions presented in the document called "*STEP\_BY\_STEP\_GUIDE.md*".<br>

Enjoy!

# References
[1]  M. A. J. Chaplain, M. Ganesh, and I. G. Graham, “Spatio-temporal pattern formation
on spherical surfaces: numerical simulation and application to solid tumour growth,”
Journal of mathematical biology, vol. 42, pp. 387–423, 2001. <br>


