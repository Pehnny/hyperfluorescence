#########################################################################################################
#
#   author(s) : Théo Piron
#   last update : 11/05/2022 (dd/mm/yyyy)
#   python version : 3.10.4
#   modules : reaseau, molecule, event, mu
#
#
#########################################################################################################
"""Module contenant les dictionnaires, fonctions et classes représentants les points et les événements.

Dictionnaries
-------------
Les valeurs sont arbitraires mais uniques pour éviter toute utilisation indésirable.
PARTICULES : {electron : 1, hole : 2, exciton : 3}
    Dictionnaire contenant les types de charges.
EVENTS : dict[str, int] = {"move" : 1, "bound" : 2, "ISC" : 3, "Forster" : 4, "decay" : 5, "unbound" : 6, "capture" : 7}
    Dictionnaire contenant les types d'événements.

Classes
-------
Point() : x, y, z
    Classe représentant un point. Les points peuvent s'additionner et se soutraire.
    Dans ce cas, le point est converti en vecteur, même dans le cas d'opération avec des nombres.
Vector(Point) : x, y, z
    Classe représentant un vecteur. Les vecteurs peuvent également se multiplier.
    La multiplication entre deux vecteurs donne le produit scalaire.
    La multiplication par un nombre donne un nouveau vecteur.

Event() : initial, final, tau, kind, particule
    Classe représentant un événement. Ceux-ci sont caractérisé par une position initiale,
    une position finale, une durée (tau), un type et une particule.
    Les opérations de comparaisons <, >, <= et >= comparent la durée des événements.
    Le but est de choisir le plus rapide.
    Les opérations de comparaion == et != comparent les événéments les autres caractéristiques des évenements.
    Le but est de déceler les actions d'une même particules ou qui causeraient une collision.
"""

from dataclasses import dataclass
from math import sqrt


EVENTS : dict[str, int] = {
    "move" : 1,
    "bound" : 2,
    "ISC" : 3,
    "ET" : 4,
    "decay" : 5,
    "capture" : 6,
}

PARTICULES : dict[str, int] = {
     "electron" : 1,
     "hole" : 2,
     "exciton" : 3
}


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
        if isinstance(other, Point) :
            return Vector(self.x + other.x, self.y + other.y, self.z + other.z)
        elif isinstance(other, (float, int)) :
            return Vector(self.x + other, self.y + other, self.z + other)
        raise TypeError(f"other must be of type Vector, float or int, got {type(other)}")
            
    def __sub__(self, other) :
        if isinstance(other, Point) :
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
        if isinstance(other, Point) :
            return self.x * other.x + self.y * other.y + self.z * other.z
        elif isinstance(other, (float, int)) :
            return Vector(self.x * other, self.y * other, self.z * other)
        raise TypeError(f"other must be of type Vector, float or int, got {type(other)}")
    
    def __rmul__(self, other) :
        return self * other
    
    def norm(self) -> float :
        return sqrt(self.x**2 + self.y**2 + self.z**2)


@dataclass(eq = False)
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
    kind : int
        Type d'événement. Les valeurs possibles sont stockées dans le EVENTS.
    particule : int = 0
        Type de particule impliquée par l'événement. Les valeurs possible sont stockées dans PARTICULES.
        La valeur par défaut peut être utilisée pour des événements spéciaux n'impliquant pas de particule.
    """

    initial : Point
    final : Point
    tau : float
    kind : int
    particule : int = 0

    def __eq__(self, other) -> bool :
        if isinstance(other, Event) :
            if self.kind == other.kind and self.particule == other.particule :
                if self.kind == EVENTS["move"] :
                    return self.initial == other.initial or self.final == other.final
                else :
                    return self.initial == other.initial and self.final == other.final
            else :
                return False
        raise TypeError(f"other must be of type event, got {type(other)}")
    
    def __ne__(self, other : object) -> bool:
        return not self == other
            
    def __lt__(self, other) -> bool :
        if isinstance(other, Event) :
            return self.tau < other.tau
        raise TypeError(f"other must be of type event, got {type(other)}")
    
    def __gt__(self, other) -> bool :
        if isinstance(other, Event) :
            return self.tau > other.tau
        raise TypeError(f"other must be of type event, got {type(other)}")

    def __le__(self, other) -> bool :
        if isinstance(other, Event) :
            return self.tau <= other.tau
        raise TypeError(f"other must be of type event, got {type(other)}")
    
    def __ge__(self, other) -> bool :
        if isinstance(other, Event) :
            return self.tau >= other.tau
        raise TypeError(f"other must be of type event, got {type(other)}")
