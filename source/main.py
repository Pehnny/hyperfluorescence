from reseau import *
from pathlib import Path



def main() -> int :
    source : Path = cwd()
    try :
        with open(source.joinpath("in.txt")) as file :
            line : str = file.readline()
    except FileNotFoundError :
        return 1
    proportions : list[float] = [float(x) for x in line.split()] 
    try :
        oled : Lattice = Lattice(proportions)
    except IndexError :
        return 2
    except ValueError :
        return 3
    value = IQE(oled)
    with open(source.joinpath("out.txt"), "w") as file :
        file.write(str(value))
    return 0

def errors(value : int) -> None :
    source : Path = cwd()
    if value != 0 :
        home : Path = source.parent
        home.joinpath("STOP").touch()
    prompt_message : str = ""
    if value == 1 :
        prompt_message += "FileNotFoundError encountered when trying to read in.txt."
    elif value == 2 :
        prompt_message += "IndexError encountered when trying to create a Lattice object."
    elif value == 3 :
        prompt_message += "ValueError encountered when trying to create a Lattice object."
    with open(source.joinpath("errors.txt"), "w") as file :
        file.write(prompt_message)

def cwd() -> Path :
    return Path(__file__).parent

def IQE(lattice : Lattice) -> float :
    recombinations : int = 10**3
    lattice.operations(recombinations)
    return 100. - lattice.get_IQE()



if __name__ == "__main__" :
    exit_code = main()
    errors(exit_code)
else :
    source : Path = cwd()
    home : Path = source.parent
    home.joinpath("STOP").touch()
    with open(source.joinpath("errors.txt"), "w") as file :
        file.write(f"Expected __name__ to be __main__, got {__name__}")
