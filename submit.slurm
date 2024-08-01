#!/bin/bash
#SBATCH --job-name=bin_inj
#SBATCH --nodes=1
#SBATCH --ntasks=96
#SBATCH --cpus-per-task=1
#SBATCH	--time=23:01:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --output=logs/R-%x.%j.out
#SBATCH --error=logs/R-%x.%j.err
#SBATCH --mail-user=wcukier.notify@gmail.com
#SBATCH --array=53-54
module purge
module load anaconda3/2022.5



python run_planets.py configs/known_systems/earth.json "no_inj" 0 $SLURM_ARRAY_TASK_ID



