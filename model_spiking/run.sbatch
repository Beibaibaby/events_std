#!/bin/bash
#SBATCH --job-name=draco            # Job name
#SBATCH --account=doiron-lab               # Account under which to submit
#SBATCH --partition=tier1q               # Partition/queue to submit to
#SBATCH --nodes=1                        # Number of nodes
#SBATCH --ntasks-per-node=1              # Number of tasks per node
#SBATCH --time=00:30:30                  # Time limit
#SBATCH --cpus-per-task=1               # Number of CPU cores per task
#SBATCH --mem=80gb                      # Memory per node
#SBATCH --output=/gpfs/data/doiron-lab/draco/out/job_output.%j.txt  # Standard output file
#SBATCH --error=/gpfs/data/doiron-lab/draco/out/job_error.%j.err    # Standard error file

# Execute the Julia script with parameters
julia run.jl --T 300 --c_noise 1.0 --sigma_noise 1.0 --ie_sign false --ee_sign false --d_ee 0.98 --f_ee 0.01 --d_ie 0.001 --f_ie 0.0 --stim_start_2 1000 --stimstr_2 -100