#!/bin/bash
#SBATCH --job-name=draco_test
#SBATCH --account=dracoxu
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:30:30
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb
#SBATCH --output=/gpfs/data/doiron-lab/draco/trash/jobname.%A_%a.out
#SBATCH --error=/gpfs/data/doiron-lab/draco/trash/jobname.%A_%a.err
#SBATCH --array=1-5%5

# Define an array of f_ie values
f_ie_values=(0.7 0.8 0.85 0.9 0.95)

# Get the f_ie value for the current task
f_ie_index=$(($SLURM_ARRAY_TASK_ID - 1))
f_ie=${f_ie_values[$f_ie_index]}

# Fixed parameters
d_ee=0.0005
f_ee=0.0
d_ie=0.15
stimstr=0.7
T=1000
stimstr_2=0
# Run the Julia script with the current parameters
julia run_4.0.jl --d_ee $d_ee --f_ee $f_ee --d_ie $d_ie --f_ie $f_ie --stimstr $stimstr --T $T --stimstr_2 $stimstr_2
