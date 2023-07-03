"""
   run.py : main gate for a multi-agent implementation of CMA-ES

   Author : Alexandre Mayer, University of Namur, Belgium
"""
from supervisor import supervisor
from worker import worker
from pathlib import Path
import argparse



def main(job_type : str, id_number : int) -> int :
    print("Python code received job type = ", job_type, " and ID = ", id_number)
    if job_type not in ["s", "w"] :
        return 1
    if not isinstance(id_number, int) :
        return 2
    if job_type == "w" :
        worker(id_number)
    elif job_type == "s" :
        supervisor(id_number)
    return 0
    
def errors(value : int) -> None :
    home : Path = cwd()
    if value != 0 :
        home.joinpath("STOP").touch()
    prompt_message : str = ""
    if value == 1 :
        prompt_message += "Unknown job type encountered."
    elif value == 2 :
        prompt_message += "Unvalid job id encountered."
    with open(home.joinpath("errors.txt"), "w") as file :
        file.write(prompt_message)
        
def cwd() -> Path :
   return Path(__file__).parent



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("job", action = "store", type = str, nargs = "?")
    parser.add_argument("id", action = "store", type = int, nargs = "?")
    args = parser.parse_args()
    job_type, id = args.job, args.id
    exit_code = main(job_type, id)
    errors(exit_code)
else :
    home = cwd()
    home.joinpath("STOP").touch()
    with open(home.joinpath("errors.txt"), "w") as file :
        file.write(f"Expected __name__ to be __main__, got {__name__}")