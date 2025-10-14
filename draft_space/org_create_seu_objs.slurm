#!/bin/bash
#SBATCH --partition="sapphire"
#SBATCH --nodes=1
#SBATCH --account="punim2251"
#SBATCH --ntasks=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=mschauhan@student.unimelb.edu.au
#SBATCH --cpus-per-task=8
#SBATCH --mem=300gb
#time in days-hours:mins:sec
#SBATCH --time=0-03:30:00

# Load modules
module load GCC/11.3.0
module OpenMPI/4.1.4
module load R/4.4.0

Rscript /data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space/ch3_qc_flames_v220_p2.R