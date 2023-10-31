# Update the .bashrc so that FEniCS only uses one core at a time
echo '# Make sure FEniCS only uses one core at a time' >> ~/.bashrc
echo 'export OMP_NUM_THREADS=1' >> ~/.bashrc 
# Reload the ./bashrc script
source ~/.bashrc
# Enter the Code folder
cd ./Code
# Plot the perturbed eigenvalues
echo "Plot perturbed eigenvalues"
python find_parameters_and_plot_perturbed_eigenvalues.py
echo "Done"
# Initiate conda
echo "Initiate conda for you shell"
source /home/johannes/anaconda3/etc/profile.d/conda.sh
conda init bash
echo "Done"
# Reload the ./bashrc script
#source ~/.bashrc
# Generate meshes using Gmsh
echo "Activate the conda-environment gmsh_latest_version"
conda activate gmsh_latest_version
echo "Done"
echo "Generate meshes with a single hole in them"
python generate_spherical_meshes_with_holes.py
echo "Done"
echo "Deactivate conda environment"
conda deactivate
echo "Done"
echo "Activate the conda environment fenicsproject"
conda activate fenicsproject
echo "Done"
echo "Convert meshes from msh format to xdmf"
python convert_meshes_from_msh_to_xdmf.py
echo "Done"
echo "Now, we will launch all simulations. Note that it will take hours to days for these simulations to finish."
python launch_simulations_Schnakenberg_sphere_with_holes.py
echo "Done"
echo "Next, we will plot the spectral decomposition of a single iteration on the mesh without a hole, i.e. the mesh corresponding to the unit sphere."
python final_concentration_decomposition.py
echo "Done"
echo "Plot the spectral decomposition as a function of the hole radius."
python plot_perturbed_eigenfunctions.py
echo "Done"
echo "Run the data analysis where we plot the quantitative properties as function of the hole radius"
python data_analysis_of_spatial_patterns.py
echo "Done"
echo "Run the data analysis of the noise"
python slope_test.py
echo "Done"
echo "Deactivate the conda environment fenicsproject"
conda deactivate
cd ..
echo "Now, we are finished!"
