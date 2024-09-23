#!/bin/bash
#SBATCH --job-name=draco_test
#SBATCH --account=dracoxu
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:30:30
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb
#SBATCH --output=/gpfs/data/doiron-lab/draco/out/demo_job.%j.out
#SBATCH --error=/gpfs/data/doiron-lab/draco/out/demo_job.%j.err

# Run the Julia script
julia run_7_t150.jl --T 2500 --dir_name_in "/gpfs/data/doiron-lab/draco/results_corr/5.4-real10" --ie_sign true --ee_sign true --add_noise false --stimstr 5.4

#julia run_7_t150.jl --T 2500 --dir_name_in "/gpfs/data/doiron-lab/draco/results_corr/0.0-real" --ie_sign true --ee_sign true --add_noise false --stimstr 0.0

#julia run_7_t150.jl --T 2500 --dir_name_in "/gpfs/data/doiron-lab/draco/results_corr/0.0-np" --ie_sign false --ee_sign false --add_noise false --stimstr 0.0