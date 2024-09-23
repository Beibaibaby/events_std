#!/bin/bash
#SBATCH --job-name=draco_test
#SBATCH --account=dracoxu
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:30
#SBATCH --cpus-per-task=1
#SBATCH --mem=90gb
#SBATCH --output=/gpfs/data/doiron-lab/draco/out_out/demo_job.%j.txt
#SBATCH --error=/gpfs/data/doiron-lab/draco/out1/demo_job.%j.err

# Run the Julia script
julia compute_corr.jl