#########################################################################################################
#
#   author(s) : Théo Piron
#   last update : 11/05/2022 (dd/mm/yyyy)
#   python version : 3.10.4
#   modules : reaseau, molecule, event, mu
#
#
#########################################################################################################

from math import sqrt

# Classe de type point qui permet de stocker des triplets de nombres correspondant aux coordonnées cartésiennes à 3 dimensions. 
# Utilisée dans le réseau et les molécules.
class point :
    def __init__(self, x : int, y : int, z : int) :
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self) -> str:
        return str((self.x, self.y, self.z))

    def __str__(self) -> str:
        return self.__repr__()

    def __eq__(self, other):
        if not isinstance(other, point):
            return NotImplemented
        else :
            cond = self.x == other.x and self.y == other.y and self.z == other.z
            return cond
    
    def __mul__(self, other) :
        if isinstance(other, point):
            return self.x * other.x + self.y * other.y + self.z * other.z
        elif isinstance(other, (float, int)):
            return point(other * self.x, other * self.y, other * self.z)
        else :
            return NotImplemented
    
    def __rmul__(self, other) :
        return self*other
            
    def __sub__(self, other) :
        if not isinstance(other, point):
            return NotImplemented
        else :
            return point(self.x - other.x, self.y - other.y, self.z - other.z)
    
    def __add__(self, other) :
        if not isinstance(other, point):
            return NotImplemented
        else :
            return point(self.x + other.x, self.y + other.y, self.z + other.z)

    def __pow__(self, other) :
        if isinstance(other, int) and other == 2 :
            return sqrt(self.x ** other + self.y ** other + self.z ** other)
        else :
            return NotImplemented
# Classe de type event qui permet au réseau d'identifier le type d'événement microscopique, le temps de réaction et les cases concernées.
class event :
    def __init__(self, initial : point, final : point, tau : float, reaction : str) :
        self.initial = initial
        self.final = final
        self.tau = tau
        self.REACTION = reaction

    def __eq__(self, other):
        if not isinstance(other, event):
            return NotImplemented
        else :
            equal = self.initial == other.initial and self.final == other.final and self.REACTION == other.REACTION
            return equal

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
        string = (str(self.initial), str(self.final), str(self.tau), self.REACTION)
        return (str(string))
