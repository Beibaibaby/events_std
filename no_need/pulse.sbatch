#!/bin/bash
#SBATCH --job-name=draco_test
#SBATCH --account=dracoxu
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=04:30:30
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb
#SBATCH --output=/gpfs/data/doiron-lab/draco/out/jobname.%A.out
#SBATCH --error=/gpfs/data/doiron-lab/draco/out/jobname.%A.err


# Fixed parameters
stimstr=0.0

# Run the Julia script with the parameters
 julia run_7_t150.jl --stimstr 5.0 --use_init_weights true --T 800 --stim_duration 30 --corr_sign false