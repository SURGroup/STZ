#!/bin/bash -l

#SBATCH
#SBATCH --job-name=Parameter_Space_Analysis
#SBATCH --time=00-10:00:00
#SBATCH --nodes=1
#SBATCH	--ntasks-per-node=1
#SBATCH --partition=shared
#SBATCH	--mail-type=end
#SBATCH --mail-user=d.alixwill@jhu.edu
#SBATCH --account=mshiel10

module load continuum_model/2.0.0
time python analyzeParameterSpace.py 
