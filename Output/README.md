# Typical output from the simulations
*Date:* 2021-09-20<br>
*Written by:* Johannes Borgqvist<br>
	When all simulations have been launched (see the file "../STEP\_BY\_STEP\_GUIDE.md" for details), the resulting datafiles are stored in subfolders in this folder. The names of these subfolders contain all the information that characterises the particular simulation that was launched. For example, the datafiles that are stored in the subfolder that is named as follows:<br>

*h\_0\_a\_0p2\_b\_1p0\_d\_20p0\_gamma\_6p873\_sigma\_0p0001\_T\_50\_ICs\_around\_steady\_states*<br>

. The datafiles stored in this folder correspond to a simulation on the hole without a hole ("h\_0"), the parameters of the corresponding simulations were (a,b,d,gamma)=(0.2,1.0,20.0,6.873), the perturbation of the steady state concentrations which determined the initial conditions were characterised by sigma=0.0001, the final time point of the simulation was chosen to T=50, and the initial conditions of these simulations were chosen so that the concentration profile of both states in each node of the mesh were set to the respective steady state concentration plus a small perturbation determined by sigma. 

Another typical file name would be:<br>

*h\_1\_r\_0p05\_a\_0p2\_b\_1p0\_d\_20p0\_gamma\_6p873\_sigma\_0p0001\_T\_50\_ICs\_around\_steady\_states*<br>


and here the only difference between the datafiles stored in this subfolder compared to the previous one is that these simulations were launched on the mesh with a single hole located at the south pole of radius r=0.05 ("h\_1\_r\_0p05"). 

Moreover, within these folders, there are subfolders for all the repeated runs that were launched on the particular mesh. In the case where we repeat each simulation 20 times, there will be 20 subfolders named as follows: "*iteration\_0*", "*iteration\_1*",...,"*iteration\_19*". Within each of these subfolders, there are 10 datafiles stored, and these files are the following:

1. *final\_timestep\_mesh.xml*: so that the script "../Code/*data\_analysis\_of\_spatial\_patterns.py*" can read the mesh for the current simulation,
2. *final\_timestep\_u.xml*: so that the script "../Code/*data\_analysis\_of\_spatial\_patterns.py*" can read the concentration profile of the active component u of the Schnakenberg model at time t=50,
3. *final\_timestep\_v.xml*: so that the script "../Code/*data\_analysis\_of\_spatial\_patterns.py*" can read the concentration profile of the inactive component v of the Schnakenberg model at time t=50,       
4. *spectral\_coefficients.csv*: the values of the coefficient in the spectral decomposition of the active component u at time t=50 which is accessed by the script "../Code/plot\_perturbed\_eigenfunctions.py",
5. *u000000.vtu*: the initial condition of the active component u at time t=0 which can be visualised in ParaView,
6. *u000001.vtu*: the concentration of the active component u at time t=50 which can be visualised in ParaView,
7. *u.pvd*: a collection of the two previous vtu files which can be visualised in ParaView,
8. *v000000.vtu*: the initial condition of the inactive component v at time t=0 which can be visualised in ParaView,
9. *v000001.vtu*: the concentration of the inactive component v at time t=50 which can be visualised in ParaView,
10. *v.pvd*: a collection of the two previous vtu files which can be visualised in ParaView.

The vtu- and pvd-files can be opened and visualised in ParaView. The xml-files can be read by FEniCS which is necessary in order to access and analyse the data  of the corresponding simulation. Lastly, the csv-file containing the spectral decomposition can be read by [Pandas](https://pandas.pydata.org/) which is used to generate the plot of the spectral decomposition of the active component u at time t=50 as a function of the hole radius.  



