"""
   worker.py : worker procedure responsible for individual fitness evaluations

   Author : Alexandre Mayer, University of Namur, Belgium
"""
from pathlib import Path
from subprocess import call

def worker(id_number : int) -> None :
    print(f"Worker with ID = {id_number}")
    home : Path = cwd()
    workdir : Path = home.joinpath(f"worker_{id_number}")
    main : Path = workdir.joinpath(f"main.py")
    call(["python3", main])
    print(f"Worker {id_number} did its job.")

def cwd() -> Path :
    return Path(__file__).parent
