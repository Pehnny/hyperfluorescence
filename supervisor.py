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
from numpy import ndarray, int64, float64



FILES = {
    "in" : "in.json",
    "out" : "out.json",
    "solver" : "solver.pkl",
    "par" : "CMAES.json",
    "history" : "history.json"
}

def supervisor(id_number : int) -> None :
    print(f"Supervisor with ID = {id_number}")

    #   Generates useful paths and gets parameters
    home : Path = cwd()
    source : Path = home.joinpath("source")
    sources : list[Path] = [source.joinpath(file) for file in source.iterdir() if file.suffix == ".py"]
    with open(home.joinpath(FILES["par"])) as file :
        parameters : dict = json.load(file)
    workers : list[Path] = [home.joinpath(f"worker_{i+1}") for i in range(parameters["population"])]
    solver_file : Path = home.joinpath(FILES["solver"])
    global_history : Path = home.joinpath(FILES["history"])

    #   Deals with the first supervisor
    if id_number == 1 :
        #   Creates solver
        solver : CMAES = CMAES(parameters["initial"], parameters["sigma"], parameters["options"])
        population : list = solver.ask(parameters["population"])
        #   Copies files from source to worker_id and create input files (first generation)
        for n, worker in enumerate(workers) :
            worker.mkdir(exist_ok = True)
            for file in sources :
                copy(file, worker)
            input_file : Path = worker.joinpath(FILES["in"])
            with open(input_file, "w") as file :
                json.dump(list(population[n]), file)
        #   Save solver state
        with open(solver_file, "wb") as file :
            file.write(solver.pickle_dumps())
        #   Creates history
        with open(global_history, "w") as file :
            json.dump(dict(), file)
        for worker in workers :
            local_history : Path = worker.joinpath(FILES["history"])
            with open(local_history, "w") as file :
                json.dump(dict(), file)
        print(f"Supervisor {id_number} did its job.") 

    #   Deals with next supervisors
    elif 2 <= id_number <= parameters["generation"] :
        #   Gets the inputs and outputs from the previous generation
        previous_population : list = []
        fitness : list = []
        for n, worker in enumerate(workers) :
            input_file : Path = worker.joinpath(FILES["in"])
            with open(input_file) as file :
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
            #   Saves local history
            local_history : Path = worker.joinpath(FILES["history"])
            population_dict = {
                f"gen_{id_number-1}" : [*previous_population[-1], fitness[-1]] 
            }
            with open(local_history) as file :
                tmp = json.load(file)
            tmp |= population_dict
            with open(local_history, "w") as file :
                json.dump(tmp, file, sort_keys = True, indent = 4)
        #   Loads and updates solver
        with open(solver_file, "rb") as file :
            solver : CMAES = pickle.load(file)
        solver.tell(previous_population, fitness)
        #   Saves history
        update : dict = {f"gen_{id_number-1}" : convert(solver.result._asdict())}
        with open(global_history) as file :
            result = json.load(file)
        result |= update
        with open(global_history, "w") as file :
            json.dump(result, file, sort_keys = True, indent = 4)
        #   Generates a new generation
        if solver.stop() :
            home.joinpath("STOP").touch()
            print(f"Solver excited early. {solver.stop()}")
            return
        population = solver.ask(parameters["population"])
        for worker, individual in zip(workers, population) :
            input_file : Path = worker.joinpath(FILES["in"])
            with open(input_file, "w") as file :
                json.dump(list(individual), file)
        #   Saves global solver
        with open(solver_file, "wb") as file :
            file.write(solver.pickle_dumps())
        print(f"Supervisor {id_number} did its job.")

    #   Deals with the last supervisor
    elif id_number > parameters["generation"] :
        home.joinpath("STOP").touch()
        #   Gets the inputs and outputs from the last generation
        previous_population : list = []
        fitness : list = []
        for n, worker in enumerate(workers) :
            input_file : Path = worker.joinpath(FILES["in"])
            with open(input_file) as file :
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
            local_history : Path = worker.joinpath(FILES["history"])
            population_dict = {
                f"gen_{id_number-1}" : [*previous_population[-1], fitness[-1]] 
            }
            with open(local_history) as file :
                tmp = json.load(file)
            tmp |= population_dict
            with open(local_history, "w") as file :
                json.dump(tmp, file, sort_keys = True, indent = 4)
        #   Loads and updates solver
        with open(solver_file, "rb") as file :
            solver : CMAES = pickle.load(file)
        solver.tell(previous_population, fitness)
        #   Saves global history
        update : dict = {f"gen_{id_number-1}" : convert(solver.result._asdict())}
        with open(global_history) as file :
            result = json.load(file)
        result |= update
        with open(global_history, "w") as file :
            json.dump(result, file, sort_keys = True, indent = 4)
        #   Cleans workers
        # for worker in workers :
        #     clean_worker(worker)
        #   Stops CMAES
        print(f"Supervisor {id_number} finished the job.")
        # print(f"All workers were cleared. Results should be stored in {global_history}")
    else :
        print(f"Unvalid supervisor ID encountered.")
        home.joinpath("STOP").touch()

def cwd() -> Path :
    return Path(__file__).parent

def convert(result : dict) -> dict :
    for key in result.keys() :
        if isinstance(result[key], ndarray) :
            result[key] = list(result[key])
        elif isinstance(result[key], int64) :
            result[key] = int(result[key])
        elif isinstance(result[key], float64) :
            result[key] = float(result[key])
    if "stop" in result :
        result["stop"] = dict(result["stop"])
    return result

def clean_worker(worker : Path) -> None :
    for path in worker.iterdir() :
        if path.is_file() :
            path.unlink()
        elif path.is_dir() :
            clean_worker(path)
