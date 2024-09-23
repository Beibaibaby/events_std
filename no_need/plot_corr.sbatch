#!/bin/bash
#SBATCH --job-name=draco_test
#SBATCH --account=dracoxu
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:30:30
#SBATCH --cpus-per-task=10
#SBATCH --mem=100gb
#SBATCH --output=/gpfs/data/doiron-lab/draco/out/a1.%A_%a.out
#SBATCH --error=/gpfs/data/doiron-lab/draco/out/a1.%A_%a.err
#SBATCH --array=1-20%20

# Use SLURM_ARRAY_JOB_ID as the identifier
ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}

export QT_QPA_PLATFORM='offscreen'

# Generate a directory name using the array job ID and task ID
dir_name_in="/gpfs/data/doiron-lab/draco/results_new/exp_${ARRAY_JOB_ID}/${SLURM_ARRAY_TASK_ID}"

event_thre=0.8

# Create the directory along with any necessary parent directories
mkdir -p "$dir_name_in"

# Run the Julia script with the current values and pass the generated directory name
julia run_7.jl --T 1000 --dir_name_in "$dir_name_in" --c_noise 0.01 --sigma_noise 1.0 --event_thre $event_thre --large_peak_mean 100 --peak_ratio 5000 --add_noise true #--f_ee 0.0
#--corr_sign false
#--f_ee 0.3 --d_ee 0.6 --f_ee 0.025 --corr_sign false
#--corr_sign false

#--f_ee 0.0 --f_ie 0.0 --d_ee 1.0 --d_ie 1.0 --f_ee 0.0

# Create a file indicating the completion of this particular task

main_dir="/gpfs/data/doiron-lab/draco/results_new/exp_${ARRAY_JOB_ID}"

julia plot_parallel.jl --main_dir "$main_dir" --event_thre $event_thre
