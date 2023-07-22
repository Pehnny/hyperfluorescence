job_type=$1                         # s=supervisor / w=worker

declare -i ID=$2                    # supervisor / worker ID

name="CMAHF"                        # nom generique pour les jobs a soumettre

supervisor_time="00:10:00"
declare -i supervisor_mem=2048

worker_time="12:00:00"
declare -i worker_mem=4096

rm -f STOP

if [ $job_type == "s" ]; then
  sbatch  --time=${supervisor_time} --mem-per-cpu=${supervisor_mem} --job-name=${name} --output=job_${ID}.out run.sh s ${ID}
  echo Job of type supervisor with ID = ${ID} was submitted
fi

if [ $job_type == "w" ]; then
  sbatch  --time=${worker_time} --mem-per-cpu=${worker_mem} --job-name=${name} run.sh w ${ID}
  echo Job of type worker with ID = ${ID} was submitted
fi
