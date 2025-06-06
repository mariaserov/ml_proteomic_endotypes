#PBS -l walltime=6:00:00
#PBS -l select=1:ncpus=1:mem=100gb
#PBS -N CGGMCombineSortK2

cd /rds/general/project/hda_24-25/live/TDS/Group09/Consensus_Pipeline_RecodedResults/CGGM_array_job_k2

module load anaconda3/personal
source activate r413

python 2CombineMat.py

