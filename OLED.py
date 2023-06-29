from reseau import *

def IQE(proportions : list[float]) -> float :
    oled : Lattice = Lattice(proportions)
    recombinations : int = 10**3
    oled.operations(recombinations)
    print(f"injections : {oled._injection}")
    print(f"recombinations : {oled._recombination}")
    print(f"emissions : {oled._emission}")
    print(f"IQE : {oled.get_IQE()}")
    return 100. - oled.get_IQE()
