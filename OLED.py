from time import time
from reseau import Lattice
from plot import plot


dimensions = (10,10,5)
proportions = (0.84,0.15,0.01)
OP = 10**2
start = time()
test = Lattice(dimensions, proportions, charges = 4)
test.operations(OP)
# for i in range(OP) :
#     electrons, holes, excitons = test.get_particules_positions()
#     plot(electrons, holes, excitons, *dimensions, str(i))
#     test.operations(1)
# electrons, holes, excitons = test.get_particules_positions()
# plot(electrons, holes, excitons, *dimensions, str(OP))
end = time()
print(f"IQE : {test.get_IQE()}")
print(f"emissions : {test._emission}")
print(f"injections : {test._injection}")
print(f"{end - start} s")
