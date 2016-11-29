#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --partition=commons
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun ~/lammps_awsemmd_20161125/bin/lmp_serial -in PROTEIN.in
