"""
   worker.py : worker procedure responsible for individual fitness evaluations

   Author : Alexandre Mayer, University of Namur, Belgium
"""
from pathlib import Path

def worker(id : int) -> None :
    print("Worker with ID = ", id)
    home : Path = cwd()
    workdir : Path = home.joinpath(f"worker_{id}")
    main : Path = workdir.joinpath(f"main.py")
    exec(main)
    print("Work done.")

def cwd() -> Path :
    return Path(__file__).parent
