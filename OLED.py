from time import time
from reseau import *
# from plot import plot
# from pathlib import Path

# home = Path.cwd()
# savepath = home.joinpath("pictures")
proportions = (0.84,0.15,0.01)
recombinations = 10**3
# OP = 10**3
start = time()
test = Lattice(proportions)
test.operations(recombinations)
# dimensions = test.get_dimensions()
# for i in range(OP) :
#     electrons, holes, excitons = test.get_particules_positions()
#     plot(electrons, holes, excitons, *dimensions, savepath.joinpath(str(i)).with_suffix(".png"))
#     test.operations(recombinations, 1)
# electrons, holes, excitons = test.get_particules_positions()
# plot(electrons, holes, excitons, *dimensions, savepath.joinpath(str(OP)).with_suffix(".png"))
end = time()
print(f"IQE : {test.get_IQE()}")
print(f"recombinations : {test._recombination}")
print(f"emissions : {test._emission}")
print(f"injections : {test._injection}")
print(f"{end - start} s")
