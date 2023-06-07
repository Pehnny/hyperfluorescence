#########################################################################################################
#
#   author(s) : Théo Piron
#   last update : 11/05/2022 (dd/mm/yyyy)
#   python version : 3.10.4
#   modules : reaseau, molecule, event, mu
#
#
#########################################################################################################
from dataclasses import dataclass

@dataclass
class point :
    """Classe représentant un point dans une grille à 3 dimensions.

    points can add and sub.

    Attributes
    ----------
    x : int
        position selon l'axe x.
    y : int
        position selon l'axe y.
    z : int
        position selon l'axe z.
    """
    x : int
    y : int
    z : int

    def __add__(self, other) :
        if isinstance(other, point):
            return point(self.x + other.x, self.y + other.y, self.z + other.z)
        else :
            raise TypeError(f"other must be of type point, got {type(other)}")
            
    def __sub__(self, other) :
        if isinstance(other, point):
            return point(self.x - other.x, self.y - other.y, self.z - other.z)
        else :
            raise TypeError(f"other must be of type point, got {type(other)}")
        

class event :
    def __init__(self, initial : point, final : point, tau : float) :
        self.initial = initial
        self.final = final
        self.tau = tau

    def __eq__(self, other):
        if isinstance(other, event):
            return self.initial == other.initial and self.final == other.final
        else :
            raise TypeError(f"other must be of type event, got {type(other)}")
            
    def __lt__(self, other) :
        if isinstance(other, event) :
            return self.tau < other.tau
        else :
            return NotImplemented
    
    def __gt__(self, other) :
        if isinstance(other, event) :
            return self.tau > other.tau
        else :
            return NotImplemented

    def __le__(self, other) :
        if isinstance(other, event) :
            return self.tau <= other.tau
        else :
            return NotImplemented
    
    def __ge__(self, other) :
        if isinstance(other, event) :
            return self.tau >= other.tau
        else :
            return NotImplemented

    def __repr__(self) -> str:
        return str(", ").join([str(self.initial), str(self.final), str(self.tau), self.REACTION])
    
a = point(2,2,2)
b = point(1,0,-1)
print(type(a))
print(type(a.x))
c = a + b
print(type(c))
print(type(c.x))