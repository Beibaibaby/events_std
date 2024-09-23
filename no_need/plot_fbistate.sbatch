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
#SBATCH --array=1-8%8

# Use SLURM_ARRAY_JOB_ID as the identifier
ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}

export QT_QPA_PLATFORM='offscreen'

# Generate a directory name using the array job ID and task ID
dir_name_in="/gpfs/data/doiron-lab/draco/results_bistate/exp_bistate_${ARRAY_JOB_ID}/${SLURM_ARRAY_TASK_ID}"

event_thre=0.8
#for periodic the parameters are Dict{String, Real}("ie_sign" => true, "add_noise" => false, "stim_duration" => 20, "T" => 1000, "scale_noise" => 1.0, "d_ee" => 0.15, "Ne" => 4000, "jie_para" => 1.4, "stim_start_2" => 400, "Nstim" => 4000, "stim_start_time" => 100, "jii_para" => -16.0, "c_noise" => 1.0, "pee" => 0.2, "pie" => 0.5, "large_peak_mean" => 300.0, "ee_sign" => true, "d_ie" => 0.25, "f_ee" => 0.3, "corr_flag" => false, "sigma_noise" => 1.0, "f_ie" => 0.3, "Ni" => 1000, "pii" => 0.5, "taui" => 10, "Ncells" => 5000, "taue" => 15, "stim_duration_2" => 200, "jee_para" => 10.0, "jei_para" => -21.599999999999998, "pei" => 0.5, "peak_ratio" => 10000.0, "stimstr_para" => 0.0, "stimstr_2" => 0.0, "K" => 800)
# Create the directory along with any necessary parent directories
mkdir -p "$dir_name_in"

# Run the Julia script with the current values and pass the generated directory name
#julia run_7_sup.jl --T 1000 --dir_name_in "$dir_name_in" --c_noise 0.01 --sigma_noise 1.0 --event_thre $event_thre --ie_sign true --ee_sign true 
julia run_7_sup.jl --T 2500 --dir_name_in "$dir_name_in" --c_noise 1.0 --sigma_noise 1.0 --event_thre $event_thre --ie_sign true --ee_sign true --d_ee 0.95 --f_ee 0.01 --d_ie 0.25 --f_ie 0.0 --stim_start_2 1000 --stimstr_2 -100
#--f_ie 0.24

# Create a file indicating the completion of this particular task

main_dir="/gpfs/data/doiron-lab/draco/results_bistate/exp_bistate_${ARRAY_JOB_ID}"

julia plot_parallel.jl --main_dir "$main_dir" --event_thre $event_thre


#julia run_7_t150.jl --T 2000 --dir_name_in "/gpfs/data/doiron-lab/draco/results_nn/crazy" --c_noise 0.01 --sigma_noise 1.0 --event_thre 0.8 --large_peak_mean 5 --peak_ratio 8000 --ie_sign false --ee_sign false

#julia plot_parallel.jl --main_dir "/gpfs/data/doiron-lab/draco/results_nn/exp_8851278" --event_thre 0.8