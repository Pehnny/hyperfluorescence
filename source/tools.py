from random import Random

class Shuffle :
    def __init__(self, array : list[int]) -> None :
        self._array : list[int] = array
        self._seed : Random = Random()
        self.length : int = len(array)

    def choice(self) -> float :
        if self.length > 0 :
            pick : int = self._seed.choice(self._array)
            self._array.remove(pick)
            self.length -= 1
            return pick
        raise IndexError("No more item in array !")

def min_0(values : list) :
    zeros = values.count(0)
    if zeros == 0 :
        return min(values)
    elif zeros < len(values) :
        subvalues = [value for value in values if value != 0]
        return min(subvalues)
    else :
        raise ValueError("Arg values contains only zeros !")
