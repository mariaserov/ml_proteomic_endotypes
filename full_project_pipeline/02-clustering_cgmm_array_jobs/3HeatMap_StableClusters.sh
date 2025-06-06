#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=1:mem=200gb
#PBS -N MakeConsPlotK2

cd /rds/general/project/hda_24-25/live/TDS/Group09/Consensus_Pipeline_RecodedResults/CGGM_array_job_k2

module load anaconda3/personal
source activate r413

python 3HeatMap_StableClusters.py

