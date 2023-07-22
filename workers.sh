#!/bin/bash
#SBATCH --job-name=CMAHF
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4096
#SBATCH --time=12:00:00
#SBATCH --array=1-50
# -------------- source exe -----------------
CURRENT=`pwd`
# ------------- copying input ---------------
# cp *.in $LOCALSCRATCH
# ----------- working directory -------------
# cd $LOCALSCRATCH
# --------- environment settings ------------
module load releases/2022a
module load Python/3.10.4-GCCcore-11.3.0   # python interpreter
# ----------------- work --------------------
echo "Job started at $(date)"
echo $CURRENT/run.py w $SLURM_ARRAY_TASK_ID
python3 $CURRENT/run.py w $SLURM_ARRAY_TASK_ID
echo "Job ended at $(date)"
# ----------- copying output ----------------
# cp *.out $CURRENT
