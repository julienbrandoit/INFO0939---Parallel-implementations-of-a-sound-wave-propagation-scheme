#!/bin/bash
#SBATCH --job-name="SERIAL output"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00
#SBATCH --output=out/serial_out1.out
#SBATCH --mem-per-cpu=4G

module load Info0939Tools

gcc SERIAL/fdtd.c -o bin/fdtd -lm -O3 
#cd ./example_inputs/slit/simple3d_4x
cd ./example_inputs/simple3d

srun ../../bin/fdtd param_3d.txt
#valgrind --tool=cachegrind ../../bin/fdtd 


#data2gmsh out_simple_4x_pz_05.dat
