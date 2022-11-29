# Remove all .mod and .o files
rm -f *.mod *.o torch*

# Remove all data
#rm -rf ./data
mkdir data
mkdir data/dynamics
#mkdir data/fields
mkdir data/densities
mkdir data/magnetization
mkdir data/lightcurves
#mkdir data/radiation
#mkdir data/orbits
#mkdir data/spectra
#mkdir data/restore
#mkdir data/particles

#rm -rf ./figs
mkdir figs

# Compile the code
module purge
#module load openmpi/openmpi-1.6.4_gcc-4.7.2_torque-4.2.3_ib
#module load hdf5/hdf5-1.8.11_szip-2.1_zlib-1.2.78_jpeglib-8d_openmpi-1.6.4_gcc-4.7.2_ib
#module load openmpi-x86_64
#module load mpi/mvapich2-x86_64
#module load mpi
#source module_load_mpi.sh
#module load libs/hdf5/1.8.19

make

# BUILD STEPS
#gfortran -I/usr/include/fgsl -g -O2 -c -o  mod_input.o mod_input.f90 
#gfortran -I/usr/include/fgsl -g -O2 -c -o  mod_circumstance.o mod_circumstance.f90 
#gfortran -I/usr/include/fgsl -g -O2 -c -o  mod_distance.o mod_distance.f90 
#gfortran -I/usr/include/fgsl -g -O2 -c -o  mod_dynamics.o mod_dynamics.f90 
#gfortran -I/usr/include/fgsl -g -O2 -c -o  mod_output.o mod_output.f90
#gfortran -I/usr/include/fgsl -g -O2 -c -o  mod_lightcurves.o mod_lightcurves.f90
#gfortran -I/usr/include/fgsl -g -O2 -c -o  main.o main.f90
#gfortran `pkg-config --cflags fgsl`  mod_input.o  mod_circumstance.o mod_distance.o  mod_dynamics.o mod_lightcurves.o mod_output.o main.o -o torch `pkg-config --libs fgsl`


# Execute the code
#qsub submit_zeltron_verus.sh
#bsub < submit_zeltron_bullet.sh
#mpirun zeltron
#mpirun -np 4 zeltron
#mpirun -n 1 -npernode 1 zeltron
#qsub submit_zeltron_psk.sh
#sbatch submit_zeltron_chuck.sh
#./torch &> ./data/torch.log
rm -f *.o *.mod *.MOD
./torch
