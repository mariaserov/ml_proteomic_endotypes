#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=4:mem=16gb
#PBS -N limma_analysis

# Change to working directory
cd /rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/limma

# Load your conda environment
source /rds/general/user/kja24/home/anaconda3/bin/activate r413

# Run the script
Rscript limma.r
