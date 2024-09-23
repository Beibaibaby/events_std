salloc --job-name=draco_interactive --account=dracoxu --partition=tier1q --nodes=1 --ntasks-per-node=1 --time=02:30:30 --cpus-per-task=1 --mem=100gb

(base) [dracoxu@cri22in002 Balanced_Sipiking]$ srun --pty $SHELL