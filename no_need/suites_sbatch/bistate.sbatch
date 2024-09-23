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
#SBATCH --array=1-6%6

# Define an array of stimstr_2 values
stimstr_2_values=(-1.7 -5.0 -10 -25 -50 0)

# Get the stimstr_2 value for the current task
stimstr_2_index=$(($SLURM_ARRAY_TASK_ID - 1))
stimstr_2=${stimstr_2_values[$stimstr_2_index]}

# Fixed parameters
d_ee=0.995
f_ee=0.0001
d_ie=0.25
f_ie=0.0

# Run the Julia script with the current parameters
julia run_4.0.jl --d_ee $d_ee --f_ee $f_ee --d_ie $d_ie --f_ie $f_ie --stimstr_2 $stimstr_2
