#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --partition=commons
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=20:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
cd /scratch/wl45/4cpv/simulation/NUMBER
srun ~/bin/lmp_serial -in 4cpv.in
