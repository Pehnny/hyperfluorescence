#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1

module load releases/2022a
module load Python/3.10.4-GCCcore-11.3.0

name="CMAES"
job_type=$1

declare -i supervisor_id=$2
declare -i nmax=101
declare -i npop=50

supervisor_time="00:10:00"
declare -i supervisor_mem=2048

worker_time="12:00:00"
declare -i worker_mem=2048

function supervisor {
    echo "I'm the supervisor #${supervisor_id} !"
    echo
    echo python3 ./run.py s ${supervisor_id}
    echo
    python3 ./run.py s ${supervisor_id}
}

function worker {
    echo "I'm the worker #${worker_id} requested by supervisor #${supervisor_id} !"
    echo
    echo python3 ./run.py w ${worker_id}
    echo
    python3 ./run.py w ${worker_id}
} 

function job {
    # Stops if "STOP" file in directory
    if [ -f "STOP" ]; then
        echo	
        echo "Stopping because a STOP signal was received !"	
        exit 0
    fi
    # Deals with a supervisor
    if [ $job_type == "s" ]; then
        supervisor
        # Stops if "STOP" file in directory
        if [ -f "STOP" ]; then
            echo	
            echo "Stopping because a STOP signal was received !"	
            exit 0
        fi
        # Stops if nmax is reached.
        if [ ${supervisor_id} -gt ${nmax} ]; then
            echo	
            echo "Stopping because the maximal number of supervisors was reached !"	
            exit 0
        fi
        # Submit a generation of workers
        out=$(sbatch --time=${worker_time} --mem-per-cpu=${worker_mem} --job-name=${name}_${supervisor_id} --array=1-${npop} $0 w ${supervisor_id})
        exit_code=$?
        if [ "$exit_code" -ne "0" ]; then
            exit 1
        fi
        # Submit next supervisor when done with workers
        tmp=${out#*Submitted batch job }
        job_id=${tmp% on cluster hercules2*}
        printf -v dependency ":"$job_id"_%i" $(seq 1 $npop)      
        out=$(sbatch --time=${supervisor_time} --mem-per-cpu=${supervisor_mem} --job-name=${name}_$(( supervisor_id + 1 )) --output=job_$(( supervisor_id + 1 )).out --dependency=afterok${dependency} $0 s $(( supervisor_id + 1 )))
        exit 0
    # Deals with a worker
    elif [ $job_type == "w" ]; then
        worker_id=${SLURM_ARRAY_TASK_ID}
        worker
        exit 0
    # Stops if neither a supervisor nor a worker was submitted
    else
        echo "Unknown type. Exciting."
        exit 1
        fi
}

# Create first supervisor or repeat itself
if [ -n "$SLURM_JOB_ID" ]; then
    job
else
    sbatch --time=${supervisor_time} --mem-per-cpu=${supervisor_mem} --job-name=${name} --output=job_1.out $0 s 1
fi
