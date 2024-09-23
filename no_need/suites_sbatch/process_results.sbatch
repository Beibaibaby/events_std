#!/bin/bash
#SBATCH --job-name=process_results
#SBATCH --account=dracoxu
#SBATCH --partition=tier1q
#SBATCH --time=01:00:00  # adjust time as needed
#SBATCH --mem=10gb       # adjust memory as needed
#SBATCH --dependency=afterok:3165411

# Read the job identifier from the file
JOB_ID=$(unique_job_id_int)

# Define the main directory where all the results are stored
main_dir="/gpfs/data/doiron-lab/draco/results_corr/exp_${JOB_ID}"

# Run the Julia script to process and plot the results
julia plot_parallel.jl --main_dir "$main_dir"
