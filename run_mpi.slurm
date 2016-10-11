#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --partition=commons
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun /home/wl45/lammps-30Jul16/src/lmp_mpi -in PROTEIN.in
