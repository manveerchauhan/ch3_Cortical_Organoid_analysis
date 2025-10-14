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
#SBATCH --mem=10gb
#time in days-hours:mins:sec
#SBATCH --time=0-03:30:00

# Load modules
module load Anaconda3/2024.02-1
echo "Modules loaded"

eval "$(conda shell.bash hook)"
conda activate FLAMES_v2.2
echo "FLAMES_v2.2 conda environment activated"

export PYTHONPATH=$PYTHONPATH:/data/gpfs/projects/punim0646/manveer/cDNA_Cupcake/sequence
export PYTHONPATH=$PYTHONPATH:/data/gpfs/projects/punim0646/manveer/cDNA_Cupcake

GTF="/data/gpfs/projects/punim2251/resoucres/gencode.v46.annotation.gtf"
genome="/data/gpfs/projects/punim2251/resoucres/true.hg38.analysisSet.fa"
isoformAnnotatedGTF="/data/gpfs/projects/punim2251/Aim1/isoquant_outputD25/OUT/OUT.transcript_models.gtf"

#filter gtf file to remove . strands 
awk '$7 == "+" || $7 == "-"' ${isoformAnnotatedGTF} > remove_unknownstrand.gtf

#run SQANTI
python3 /data/gpfs/projects/punim0646/manveer/SQANTI3-5.1.2/sqanti3_qc.py remove_unknownstrand.gtf ${GTF} ${genome}