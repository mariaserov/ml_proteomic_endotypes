#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=1:mem=50gb
#PBS -N dict

cd /rds/general/project/hda_24-25/live/TDS/Group09/Data

module load anaconda3/personal
source activate r413

Rscript match_filter_hypertension.R 