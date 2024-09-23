#!/bin/bash
#SBATCH --job-name=draco_test
#SBATCH --account=dracoxu
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:30:30
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb
#SBATCH --output=/gpfs/data/doiron-lab/draco/trash/jobname.%A_%a.out
#SBATCH --error=/gpfs/data/doiron-lab/draco/trash/jobname.%A_%a.err
#SBATCH --array=1-5%5

# Define arrays of values
jie_values=(1.0 1.1 1.5 0.5 2.0)
stimstr_values=(0.0)

# Total combinations of jie and stimstr
total_combinations=$((${#jie_values[@]} * ${#stimstr_values[@]}))

# Get indices for jie and stimstr
array_index=$(($SLURM_ARRAY_TASK_ID - 1))
jie_index=$(($array_index / ${#stimstr_values[@]}))
stimstr_index=$(($array_index % ${#stimstr_values[@]}))

# Get the jie and stimstr values for the current task
jie=${jie_values[$jie_index]}
stimstr=${stimstr_values[$stimstr_index]}

# Fixed parameters
d_ee=0.3
f_ee=0.4
d_ie=0.23
f_ie=0.3

# Run the Julia script with the current parameters
julia run_4.0.jl --d_ee $d_ee --f_ee $f_ee --d_ie $d_ie --f_ie $f_ie --stimstr $stimstr --jie $jie
