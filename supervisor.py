"""
   supervisor.py : main procedure for a multi-agent implementation of CMA-ES

   Author : Alexandre Mayer, University of Namur, Belgium

   Reference for CMA-ES : https://blog.otoro.net/2017/10/29/visual-evolution-strategies   
"""
import pickle
from cma import CMAEvolutionStrategy as CMAES
from pathlib import Path
from shutil import copy
import json



FILES = {
    "in" : "in.json",
    "out" : "out.json",
    "solver" : "solver.pkl",
    "par" : "CMAES.json",
    "history" : "history.txt"
}

def supervisor(id_number : int) -> None :
    print(f"Supervisor with ID = {id_number}")

    # Generates useful paths and gets parameters
    home : Path = cwd()
    source : Path = home.joinpath("source")
    sources : list[Path] = [source.joinpath(file) for file in source.iterdir() if file.suffix == ".py"]
    with open(home.joinpath(FILES["par"])) as file :
        parameters : dict = json.load(file)
    workers : list[Path] = [home.joinpath(f"worker_{i+1}") for i in range(parameters["population"])]
    solver_file : Path = home.joinpath(FILES["solver"])
    history : Path = home.joinpath(FILES["history"])

    # Deals with the first supervisor
    if id_number == 1 :
        # Creates solver
        solver : CMAES = CMAES(parameters["initial"], parameters["sigma"], parameters["options"])
        population : list = solver.ask(parameters["population"])
        # Copies files from source to worker_id and create input files (first generation)
        for n, worker in enumerate(workers) :
            worker.mkdir(exist_ok = True)
            for file in sources :
                copy(file, worker)
            input : Path = worker.joinpath(FILES["in"])
            with open(input, "w") as file :
                json.dump(list(population[n]), file)
        # Save solver state
        with open(solver_file, "wb") as file :
            file.write(solver.pickle_dumps())
        # Creates history
        try :
            with open(history, "x") as file :
                ...
        except FileExistsError :
            with open(history, "w") as file :
                file.write("")
        print(f"Supervisor {id_number} did its job.") 
    # Deals with next supervisors
    elif 2 <= id_number <= parameters["generation"] :
        # Gets the inputs and outputs from the previous generation
        previous_population : list = []
        fitness : list = []
        for n, worker in enumerate(workers) :
            input : Path = worker.joinpath(FILES["in"])
            with open(input) as file :
                previous_population.append(json.load(file))
            output : Path = worker.joinpath(FILES["out"])
            try :
                with open(output) as file :
                    fitness.append(json.load(file))
            except FileNotFoundError :
                home.joinpath("STOP").touch()
                with open(home.joinpath("errors.txt"), "w") as file :
                    file.write(f"Supervisor {id_number} didn't find the output of worker {n}")
                return
        # Loads and updates solver
        with open(solver_file, "rb") as file :
            solver : CMAES = pickle.load(file)
        solver.tell(previous_population, fitness)
        # Clears workers output
        for worker in workers :
            worker.joinpath(FILES["out"]).unlink()
        # Saves history
        xbest : list = list(solver.result.xbest)
        fbest : float = float(solver.result.fbest)
        values : list = [str(x) for x in xbest] + [str(fbest), "\n"]
        print(values)
        with open(history, "a") as file :
            file.write(str("\t").join(values))
        # Generates a new generation
        if solver.stop() :
            home.joinpath("STOP").touch()
            print(f"Solver excited early.")
            return
        population : list = solver.ask(parameters["population"])
        for n, worker in enumerate(workers) :
            input : Path = worker.joinpath(FILES["in"])
            with open(input, "w") as file :
                json.dump(list(population[n]), file)
        # Saves solver
        with open(solver_file, "wb") as file :
            file.write(solver.pickle_dumps())
        print(f"Supervisor {id_number} did its job.")
    # Deals with the last supervisor
    elif id_number > parameters["generation"] :
         # Gets the inputs and outputs from the last generation
        previous_population : list = []
        fitness : list = []
        for n, worker in enumerate(workers) :
            input : Path = worker.joinpath(FILES["in"])
            with open(input) as file :
                previous_population.append(json.load(file))
            output : Path = worker.joinpath(FILES["out"])
            try :
                with open(output) as file :
                    fitness.append(json.load(file))
            except FileNotFoundError :
                home.joinpath("STOP").touch()
                with open(home.joinpath("errors.txt"), "w") as file :
                    file.write(f"Supervisor {id_number} didn't find the output of worker {n}")
                return
        # Loads and updates solver
        with open(solver_file, "rb") as file :
            solver : CMAES = pickle.load(file)
        solver.tell(previous_population, fitness)
        # Saves history
        xbest : list = list(solver.result.xbest)
        fbest : float = float(solver.result.fbest)
        values : list = [str(x) for x in xbest] + [str(fbest), "\n"]
        with open(history, "a") as file :
            file.write(str("\t").join(values))
        # Cleans workers
        for worker in workers :
            for file in worker.iterdir() :
                file.unlink()
            worker.rmdir()
        # Stops CMAES
        home.joinpath("STOP").touch()
        print(f"Supervisor {id_number} finished the job.")
        print(f"All workers were cleared. Results should be stored in {history}")
    else :
        print(f"Unvalid supervisor ID encountered.")
        home.joinpath("STOP").touch()

def cwd() -> Path :
    return Path(__file__).parent