#!/bin/bash
#SBATCH --job-name=bin_inj
#SBATCH --nodes=1
#SBATCH --ntasks=96
#SBATCH --cpus-per-task=1
#SBATCH	--time=23:59:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --output=logs/R-%x.%j.out
#SBATCH --error=logs/R-%x.%j.err
#SBATCH --mail-user=wcukier.notify@gmail.com
#SBATCH --array=0-3
module purge
module load anaconda3/2022.5



python run_planets.py configs/known_systems/earth.json "detailed_bin" 2 $SLURM_ARRAY_TASK_ID



