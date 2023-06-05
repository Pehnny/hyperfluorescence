import random
from time import time_ns
from numpy.random import default_rng

length = 10**6

start = time_ns()
rn_numbers = [random.randint(0,99) for i in range(length)]
stop = time_ns()
print(f"Temps écoulé pour random : {(stop - start)/length} ns")

seed = default_rng()
start = time_ns()
np_rn_numbers = [seed.integers(0,100) for i in range(length)]
stop = time_ns()

print(f"Temps écoulé pour defaut_rng : {(stop - start)/length} ns")