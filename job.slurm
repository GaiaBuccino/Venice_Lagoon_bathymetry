#! /bin/bash
  
#SBATCH --job-name=Read_bathymetry
#SBATCH -N1 -n 1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=350gb
#SBATCH --time=00:59:00
#SBATCH --account=OGS23_PRACE_IT
#SBATCH --partition=g100_usr_prod
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gbuccino@ogs.it


date

python ./read_bathy.py

date
