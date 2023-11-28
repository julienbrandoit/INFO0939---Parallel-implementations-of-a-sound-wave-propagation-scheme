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

param_files=("convergence4.txt" "convergence5.txt" "convergence6.txt" "convergence7.txt" "convergence8.txt" "convergence9.txt") # Ajoutez vos fichiers ici
# Loop over each parameter file
for param_file in "${param_files[@]}"
do
    # Extract the number from the file name
    number=$(echo $param_file | grep -o -E '[0-9]+')

    # Set output file and output directory names using the extracted number
    output_file="fdtd_serial_${number}.out"
    output_dir="outputs_${number}"

    # Set the SLURM output file name
    # Note: Ensure that SLURM accepts this dynamic naming within the script
    # If not, you might need to create separate submission scripts for each job
    #SBATCH --output=${output_file}

    # Run the program
    cd example_inputs/simple3d
    ../../fdtd $param_file

    # Save the output
    mkdir -p ${output_dir}
    # Move output files to the output directory
    mv out_*${number}.dat ${output_dir}/
    
done