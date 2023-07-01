"""
   supervisor.py : main procedure for a multi-agent implementation of CMA-ES

   Author : Alexandre Mayer, University of Namur, Belgium

   Reference for CMA-ES : https://blog.otoro.net/2017/10/29/visual-evolution-strategies   
"""
from Utilitaries import create_directories, prepare_workers, write_inputs, get_inputs, get_outputs, clean_workers, update_sequence
import matplotlib.pyplot as plt
import numpy as np
import pickle
from cma import CMAEvolutionStrategy as CMAES
from pathlib import Path
import json



def supervisor(id : int):
    print("Supervisor with ID = ", id)

    home : Path = cwd()
    with open(home.joinpath("CMAES.json")) as file :
        content = file.read()
    parameters = json.loads(content)
    if id == 1 :
        for i in range(parameters["population"]) :
            home.joinpath(f"worker_{i}").mkdir(exist_ok = True)


    
    # Read technical parameters from CMAES.par
    # ++++++++++++++++++++++++++++++++++++++++

    save_history = True

    # Initialization
    # ++++++++++++++
    # - First supervisor : initialize directories & start CMA-ES
    # ----------------------------------------------------------
    if ID==1:
        print("Creating directories ...\n")

        create_directories (npop)  # create directories

        prepare_workers(npop)      # copy the files needed to compute the fitness

        update_sequence(npop,ID)   # update sequence of workers / supervisors

        solver = CMAES(n, popsize=npop, weight_decay=0., sigma_init = 5., sopt=sopt)
                                # Initialize CMA-ES solver
        print()

    # - Next supervisors : reload previous work
    # -----------------------------------------
    if ID>1:
        points = get_inputs (npop, n)

        list_fitness, ready = get_outputs (npop)

        if not ready:
            print("\nWorker calculations are not yet finished ; try again later !")
            file=open('STOP','w') ; file.close()
            return
        else:
            print("Worker calculations were loaded successfully.\n")

        if ID<ngen+1:
            update_sequence(npop,ID)   # update sequence of workers / supervisors

        if save_history:
            if ID==2:
                History = []
        else:
            History = pickle.load(open('saved-history.pkl','rb'))

        History.append([points,list_fitness])

        pickle.dump(History,open('saved-history.pkl','wb'))

    # CMA-ES must be restarted in its previous state

        print("Restarting CMA-ES ...\n")

        solver = pickle.load(open('saved-solver-object.pkl','rb')) ; print("CMA-ES restarted !\n")

    # Tell CMA-ES the fitness values

        solver.tell(list_fitness)  # tell solver the fitness values

    # Keep track of the best-so-far solution and save history (starting from the second supervisor)
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if ID>=2:
        if ID==2:
        xbest = points[0] ; fbest = list_fitness[0] ; List_fbest = []
        else:
        file=open('CMAES.out','r')
        xbest = np.zeros(n)
        for i in range(n):
            line=file.readline() ; xbest[i] = float(line)
        line=file.readline() ; fbest = np.array(float(line))
        file.close()

        List_fbest = pickle.load(open('saved-list-fbest.pkl','rb'))

        if sopt==1:
        for i in range(npop): 
            if list_fitness[i] > fbest:
            xbest = points[i] ; fbest = list_fitness[i]
        else:  
        for i in range(npop): 
            if list_fitness[i] < fbest:
            xbest = points[i] ; fbest = list_fitness[i]

        print("Best fitness so far : %f\n" % fbest)

        file=open('CMAES.out','w')
        for i in range(n):
        file.write(str(xbest[i])+"\n")
        file.write(str(fbest)+"\n")
        file.close()

        List_fbest.append(fbest)

        pickle.dump(List_fbest,open('saved-list-fbest.pkl','wb'))

        plt.plot(List_fbest,linewidth=3)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.xlabel("#generations",size=14)
        plt.ylabel("fitness",size=14)
        plt.title("Best fitness as a function of the number of generations",size=14)
        plt.savefig("Fitness.png")
        plt.close()

    # Ask CMA-ES the next points to evaluate
    # ++++++++++++++++++++++++++++++++++++++
    if ID <= ngen:
        print("Writing new input files ...\n")

        points = np.array(solver.ask())            # ask solver for points to evaluate

        write_inputs(points)  

    # Save CMA-ES in its current state

        pickle.dump(solver,open('saved-solver-object.pkl','wb')) ; print("CMA-ES saved\n")

    # Finish after ngen generations
    # +++++++++++++++++++++++++++++
    if ID > ngen:
        print("This is all !\n")
    
        print("Local optimum discovered by the solver : ", xbest, "\n")
        print("Fitness at this local optimum : ", fbest, "\n")
    
        file=open('STOP','w') ; file.close()

        clean_workers(npop)
        
    print("Bye, now !")

def cwd() -> Path :
    return Path(__file__).parent