#!/bin/bash
#SBATCH --job-name="fdtd serial"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00

# Compile the program
module load Info0939Tools
module load GCC
gcc -O3 -o fdtd fdtd.c -lm

# List of parameter files

param_files=("convergence3.txt" "convergence4.txt" "convergence5.txt" "convergence6.txt" "convergence7.txt") # Ajoutez vos fichiers ici
# Loop over each parameter file
for param_file in "${param_files[@]}"
do
    # Extract the number from the file name
    number=$(echo $param_file | grep -o -E '[0-9]+')

    # Set output file and output directory names using the extracted number
    output_file="fdtd_serial_${number}.out"
    output_dir="outputs_${number}"

    # Run the program
    cd example_inputs/simple3d
    ../../fdtd $param_file

    # Save the output
    mkdir -p ${output_dir}
    # Move output files to the output directory
    mv out_*${number}.dat ${output_dir}/
    
done