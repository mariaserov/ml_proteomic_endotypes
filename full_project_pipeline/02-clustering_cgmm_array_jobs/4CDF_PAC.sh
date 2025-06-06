#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=1:mem=50gb
#PBS -N CDFk2

cd /rds/general/project/hda_24-25/live/TDS/Group09/Consensus_Pipeline_RecodedResults/CGGM_array_job_k2

module load anaconda3/personal
source activate r413

python 4CDF_PAC.py
