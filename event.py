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
    """

    ...


@dataclass
class Isc(Event) :
    """Dataclasse représentant un événement de type conversion intersystème au sein du réseau.

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

    ...


@dataclass
class Decay(Event) :
    """Dataclasse représentant un événement de type décroissance au sein du réseau.

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

    ...


