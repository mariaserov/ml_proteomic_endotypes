#PBS -l walltime=6:00:00
#PBS -l select=1:ncpus=1:mem=20gb
#PBS -N PAC

cd /rds/general/project/hda_24-25/live/TDS/Group09/Consensus_Pipeline_RecodedResults

module load anaconda3/personal
source activate r413

python 8PAC.py
