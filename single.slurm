#!/bin/bash
#SBATCH --job-name=earth
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH	--time=03:59:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --output=logs/R-%x.%j.out
#SBATCH --error=logs/R-%x.%j.err
#SBATCH --mail-user=wcukier.notify@gmail.com

module purge
module load anaconda3/2022.5

## mass_2=$(python3 -c "print(1/30 * $SLURM_ARRAY_TASK_ID)")


python main.py configs/known_systems/earth_inc.json "earth_inc" 1

###python main.py configs/known_systems/earth.json q-d_1au_1me m2 $mass_2