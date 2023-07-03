from reseau import *
from pathlib import Path
import json


FILES = {
    "in" : "in.json",
    "out" : "out.json"
}

def main() -> int :
    source : Path = cwd()
    try :
        with open(source.joinpath(FILES["in"])) as file :
            input : list = json.load(file)
        proportions = (1. - sum(input), *input)
    except FileNotFoundError :
        return 1
    oled : Lattice = Lattice(proportions)
    value = IQE(oled)
    with open(source.joinpath(FILES["out"]), "w") as file :
        json.dump(value, file)
    return 0

def errors(value : int) -> None :
    source : Path = cwd()
    if value != 0 :
        home : Path = source.parent
        home.joinpath("STOP").touch()
    prompt_message : str = ""
    if value == 1 :
        prompt_message += "FileNotFoundError encountered when trying to read in.txt."
    with open(source.joinpath("errors.txt"), "w") as file :
        file.write(prompt_message)

def cwd() -> Path :
    return Path(__file__).parent

def IQE(lattice : Lattice) -> float :
    recombinations : int = 10**0
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
