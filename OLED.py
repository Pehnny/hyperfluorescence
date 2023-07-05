from reseau import *

def IQE(proportions : list[float], dimension : tuple[int, int, int], charges : int) -> float :
    oled : Lattice = Lattice(proportions, dimension=dimension, charges=charges)
    recombinations : int = 10**3
    oled.operations(recombinations)
    print(f"dimensions : {dimension}")
    print(f"charges : {charges}")
    print(f"injections : {oled._injection}")
    print(f"recombinations : {oled._recombination}")
    print(f"emissions : {oled._emission}")
    print(f"IQE : {oled.get_IQE()}")
    return 100. - oled.get_IQE()

proportions = (0.84,0.15,0.01)
dimension = (20,20,20)
charges = 4
IQE(proportions, dimension, charges)
