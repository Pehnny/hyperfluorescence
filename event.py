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
class Point :
    """Dataclasse représentant un point dans une grille à 3 dimensions.

    Les points peuvent s'additionner et se soustraire

    Attributes
    ----------
    x : float
        Position selon l'axe x.
    y : float
        Position selon l'axe y.
    z : float
        Position selon l'axe z.
    """
    x : int
    y : int
    z : int

    def __add__(self, other) :
        if isinstance(other, (Point, Vector)) :
            return Vector(self.x + other.x, self.y + other.y, self.z + other.z)
        elif isinstance(other, (float, int)) :
            return Vector(self.x + other, self.y + other, self.z + other)
        raise TypeError(f"other must be of type Vector, float or int, got {type(other)}")
            
    def __sub__(self, other) :
        if isinstance(other, (Point, Vector)) :
            return Vector(self.x - other.x, self.y - other.y, self.z - other.z)
        elif isinstance(other, (float, int)) :
            return Vector(self.x - other, self.y - other, self.z - other)
        raise TypeError(f"other must be of type Vector, float or int, got {type(other)}")


@dataclass
class Vector(Point) :
    x : float
    y : float
    z : float
    
    def __mul__(self, other) :
        if isinstance(other, (Vector, Point)) :
            return self.x * other.x + self.y * other.y + self.z * other.z
        elif isinstance(other, (float, int)) :
            return Vector(self.x * other, self.y * other, self.z * other)
        raise TypeError(f"other must be of type Vector, float or int, got {type(other)}")
    
    def __rmul__(self, other) :
        return self * other


@dataclass
class Event :
    """Dataclasse représentant un événement au sein du réseau.

    Opération de comparaison implémentée

    Attributes
    ----------
    initial : Point
        Position de départ de l'événement.
    final : Point
        Position d'arrivée de l'événement.
    tau : float
        Durée de l'événement.
    """
    initial : Point
    final : Point
    tau : float
            
    def __lt__(self, other) :
        if isinstance(other, Event) :
            return self.tau < other.tau
        raise TypeError(f"other must be of type event, got {type(other)}")
    
    def __gt__(self, other) :
        if isinstance(other, Event) :
            return self.tau > other.tau
        raise TypeError(f"other must be of type event, got {type(other)}")

    def __le__(self, other) :
        if isinstance(other, Event) :
            return self.tau <= other.tau
        raise TypeError(f"other must be of type event, got {type(other)}")
    
    def __ge__(self, other) :
        if isinstance(other, Event) :
            return self.tau >= other.tau
        raise TypeError(f"other must be of type event, got {type(other)}")


class Particules :
    """Classe contenant les types de particules pour les événements de type Move.

    Attributes
    ----------
    electron : str = "electron"

    hole : str = "hole"

    exciton : str = "exciton"    
    """

    electron : str = "electron"
    hole : str = "hole"
    exciton : str = "exciton"


@dataclass
class Move(Event) :
    """Dataclasse représentant un événement de type déplacement.

    Opération de comparaison implémentée.

    Attributes
    ----------
    initial : Point
        Position de départ de l'événement.
    final : Point
        Position d'arrivée de l'événement.
    tau : float
        Durée de l'événement.
    particule : str
        Type de particule sujette au déplacement. Les types autorisés sont "electron", "hole" et "exciton". \n
        Ces valeurs sont stockées comme attributs de classe dans la classe Particules.
    """

    particule : str


class Spins :
    """Classe contenant les sens de conversion intersystème pour les événements de type Isc.

    Attributes
    ----------
    direct : str = "direct"
    reverser : str = "reverse"
    """

    direct : str = "ISC"
    reverse : str = "RISC"


@dataclass
class Isc(Event) :
    """Dataclasse représentant un événement de type conversion intersystème.

    Opération de comparaison implémentée.

    Attributes
    ----------
    initial : Point
        Position de départ de l'événement.
    final : Point
        Position d'arrivée de l'événement.
    tau : float
        Durée de l'événement.
    conversion : str
        Sens de la conversion intersystème. Peut être "direct" (triplet-singlet) ou "reverse" (singlet-triplet). \n
        Les sens de réaction sont stockée dans la classe Spins.
    """

    conversion : str


class Radiation :
    """Classe contenant les types de décroissances pour les événements de type Decay.

    Attributes
    ----------
    non_radiative : str = "NR"

    fluorescent : str = "Fluo"
    """

    non_radiative : str = "NR"
    fluorescent : str = "Fluo"


@dataclass
class Decay(Event) :
    """Dataclasse représentant un événement de type décroissance.

    Opération de comparaison implémentée.

    Attributes
    ----------
    initial : Point
        Position de départ de l'événement.
    final : Point
        Position d'arrivée de l'événement.
    tau : float
        Durée de l'événement.
    radiative : str
        Type de radiation. Les types de radiations sont "NR" (non radiation), sans émission, \n
        et "Fluo" (fluorescence), avec émission fluorescente. Les valeurs possibles sont stockées \n
        dans la classe Radiation.
    """

    radiative : str


@dataclass
class Capture(Event) :
    """Dataclasse représentant un événement de type capture électronique.

    Opération de comparaison implémentée.

    Attributes
    ----------
    initial : Point
        Position de départ de l'événement.
    final : Point
        Position d'arrivée de l'événement.
    tau : float
        Durée de l'événement.
    particule : str
        Type de particule sujette au déplacement. Les types autorisés sont "electron", "hole" et "exciton". \n
        Ces valeurs sont stockées comme attributs de classe dans la classe Particules.
    """

    particule : str


@dataclass
class Unbound(Event) :
    """Classe représentant un événement de type séparation des charges d'un exciton."""
    ...
