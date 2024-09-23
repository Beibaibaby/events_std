#!/bin/bash
#SBATCH --job-name=draco_test
#SBATCH --account=dracoxu
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:30:30
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb
#SBATCH --output=/gpfs/data/doiron-lab/draco/trash_suites/b.%A_%b.out
#SBATCH --error=/gpfs/data/doiron-lab/draco/trash_suites/b.%A_%b.err
# Update array range based on total combinations
#SBATCH --array=1-50%25

export QT_QPA_PLATFORM='offscreen'

# Define arrays of values
stimstr_values=(10 20)
f_ees=(0.85)
f_ies=(0.0)

# Total combinations of stimstr, f_ee, and f_ie
total_combinations=$((${#stimstr_values[@]} * ${#f_ees[@]} * ${#f_ies[@]}))

# Get indices for stimstr, f_ee, and f_ie
array_index=$(($SLURM_ARRAY_TASK_ID - 1))

# Calculate indices for each parameter
stimstr_index=$(($array_index / (${#f_ees[@]} * ${#f_ies[@]})))
f_ee_index=$(($array_index / ${#f_ies[@]} % ${#f_ees[@]}))
f_ie_index=$(($array_index % ${#f_ies[@]}))

# Get the values for the current task
stimstr=${stimstr_values[$stimstr_index]}
f_ee=${f_ees[$f_ee_index]}
f_ie=${f_ies[$f_ie_index]}

# Fixed parameters
d_ee=0.24
d_ie=0.24

# Run the Julia script with the current parameters
julia run_6.jl --T 1000 --d_ee $d_ee --f_ee $f_ee --d_ie $d_ie --f_ie $f_ie --stimstr $stimstr --corr_sign false --add_noise false
