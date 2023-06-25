from time import time
from reseau import *
# from plot import plot

proportions = (0.84,0.15,0.01)
recombinations = 10**5
start = time()
test = Lattice(proportions)
test.operations(recombinations)
# dimensions = test.get_dimensions()
# for i in range(OP) :
#     electrons, holes, excitons = test.get_particules_positions()
#     plot(electrons, holes, excitons, *dimensions, str(i))
#     test.operations(1)
# electrons, holes, excitons = test.get_particules_positions()
# plot(electrons, holes, excitons, *dimensions, str(OP))
end = time()
print(f"IQE : {test.get_IQE()}")
print(f"recombinations : {test._recombination}")
print(f"emissions : {test._emission}")
print(f"injections : {test._injection}")
print(f"{end - start} s")
