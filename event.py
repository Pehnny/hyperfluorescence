#########################################################################################################
#
#   author(s) : Théo Piron
#   last update : 11/05/2022 (dd/mm/yyyy)
#   python version : 3.10.4
#   modules : reaseau, molecule, event, mu
#
#
#########################################################################################################
class point :
    def __init__(self, x : int, y : int, z : int) :
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self) -> str:
        return str((self.x, self.y, self.z))

    def __eq__(self, other):
        if isinstance(other, point):
            return self.x == other.x and self.y == other.y and self.z == other.z
        else :
            raise TypeError(f"other must be of type point, got {type(other)}")
            
    def __sub__(self, other) :
        if isinstance(other, point):
            return point(self.x - other.x, self.y - other.y, self.z - other.z)
        else :
            raise TypeError(f"other must be of type point, got {type(other)}")
    
    def __add__(self, other) :
        if isinstance(other, point):
            return point(self.x + other.x, self.y + other.y, self.z + other.z)
        else :
            raise TypeError(f"other must be of type point, got {type(other)}")
        

# Classe de type event qui permet au réseau d'identifier le type d'événement microscopique, le temps de réaction et les cases concernées.
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