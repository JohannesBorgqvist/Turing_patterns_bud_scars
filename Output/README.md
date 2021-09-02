# The output folder
Here, all the pvd files from the FEMFD-simulations of the Schnackenberg model on the sphere with holes is stored. However, these files are way too large in terms of memory to push them to github, so therefore these are added to the gitignore file. When the script is executed all sub folder in the folder "./Output" will be named based on the parameters used in that particular simulation. The sub folders are named as follows:<br>

h\_*numHoles*\_a\_*a_val*\_b\_*b_val*\_d\_*d_val*\_gamma\_*gamma_val*\_sigma\_*sigma_val*\_T\_*T_val*\_laa\_*laa_val*\_lab\_*lab_val*\_cg_*cg\_str*/<br>
Here, follows the explanation for each part of this rather long folder name:
* "h\_*numHoles*": h indicates the number of holes on the sphere where *numHoles* is an integer,
* "a\_*a_val*": a is one of the two rate parameters in the Schnackenberg model where *a_val* is a string of the numerical value rounded to two digits where the point is replaced by the character "p" (e.g. if a=0.25 then a_val="0p25"),
* "b\_*b_val*": b is one of the two rate parameters in the Schnackenberg model where *b_val* is a string of the numerical value rounded to two digits where the point is replaced by the character "p" (e.g. if b=6.33 then b_val="6p33"),
* "d\_*d_val*": d is the relative diffusion in the Schnackenberg model where *d_val* is a string of the numerical value rounded to two digits where the point is replaced by the character "p" (see the previous points for an explanation),
* "gamma\_*gamma_val*": gamma is the relative diffusion in the Schnackenberg model where *gamma_val* is a string of the numerical value rounded to two digits where the point is replaced by the character "p",
* "sigma\_*sigma_val*": sigma is the variance of the normally distributed perturbation that is added to the initial conditions where *sigma_val* is a string of the numerical value rounded to two digits where the point is replaced by the character "p",
* "T\_*T_val*": T is the time which determines the termination criteria of the FD time stepping scheme where *T_val* is a string of the numerical value rounded to two digits where the point is replaced by the character "p",
* "laa\_*laa_val*": determines the factor which the parameter a is either locally enhanced or diminished in the region adjacent to the hole (e.g. if laa=1 then we have no local activation of a, if laa=2 then a is enhanced in the adjacent region and if laa=0.5 then it is weakened in the adjacent region),
* "lab\_*lab_val*": the same as the previous string but with respect to the local activation of the rate parameter b in the Schnackenberg model,
* "cg\_*cg\_str*": where cg corresponds to whether the cell had cell growth or not and the string *cg\_str* either takes the value "yes" or "no".

Thus, the names of each sub folder contains all the information needed in order to grasp the conditions for each simulation. In each of these folders there are two main files called "u.pvd" and "v.pvd" which captures the spatiotemporal evolution of the Schnackenberg RD-model and it is easily visualised with the software "*ParaView*" which is described further in the main "README.md" file. 
