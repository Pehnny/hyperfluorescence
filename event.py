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
    x : int
        Position selon l'axe x.
    y : int
        Position selon l'axe y.
    z : int
        Position selon l'axe z.
    """
    x : int
    y : int
    z : int

    def __add__(self, other) :
        if isinstance(other, Point):
            return Point(self.x + other.x, self.y + other.y, self.z + other.z)
        raise TypeError(f"other must be of type point, got {type(other)}")
            
    def __sub__(self, other) :
        if isinstance(other, Point):
            return Point(self.x - other.x, self.y - other.y, self.z - other.z)
        raise TypeError(f"other must be of type point, got {type(other)}")


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
    """Dataclasse représentant un événement de type déplacement au sein du réseau.

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
    """Dataclasse représentant un événement de type conversion intersystème au sein du réseau.

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
    """Dataclasse représentant un événement de type décroissance au sein du réseau.

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
