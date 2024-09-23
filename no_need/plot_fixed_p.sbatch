#!/bin/bash
#SBATCH --job-name=draco_test
#SBATCH --account=dracoxu
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=06:30:30
#SBATCH --cpus-per-task=10
#SBATCH --mem=100gb
#SBATCH --output=/gpfs/data/doiron-lab/draco/out/a1.%A_%a.out
#SBATCH --error=/gpfs/data/doiron-lab/draco/out/a1.%A_%a.err
#SBATCH --array=1-150%50

# Use SLURM_ARRAY_JOB_ID as the identifier
ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}

export QT_QPA_PLATFORM='offscreen'

# Generate a directory name using the array job ID and task ID
dir_name_in="/gpfs/data/doiron-lab/draco/results_nn/exp_fixed_${ARRAY_JOB_ID}/${SLURM_ARRAY_TASK_ID}"

event_thre=0.8

# Create the directory along with any necessary parent directories
mkdir -p "$dir_name_in"

# Run the Julia script with the current values and pass the generated directory name
julia run_7_t150.jl --T 2500 --dir_name_in "$dir_name_in" --c_noise 0.01 --sigma_noise 1.0 --event_thre $event_thre --peak_ratio 18000 --ie_sign false --ee_sign false

# Create a file indicating the completion of this particular task

main_dir="/gpfs/data/doiron-lab/draco/results_nn/exp_fixed_${ARRAY_JOB_ID}"

julia plot_parallel.jl --main_dir "$main_dir" --event_thre $event_thre


#julia run_7_t150.jl --T 2000 --dir_name_in "/gpfs/data/doiron-lab/draco/results_nn/crazy" --c_noise 0.01 --sigma_noise 1.0 --event_thre 0.8 --large_peak_mean 5 --peak_ratio 8000 --ie_sign false --ee_sign false

