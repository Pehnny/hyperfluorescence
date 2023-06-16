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

Dictionnaires
-------------
Les valeurs sont arbitraires mais uniques pour éviter toute utilisation indésirable.
PARTICULES : {electron : 0, hole : 1, exciton : 2}
    Dictionnaire contenant les types de charges. Les valeurs sont arbitraires.
CROSSING : {ISC : 3, RISC : 4}
    Dictionnaire contenant les sens de conversion intersystème. Les valeurs sont arbitraires.
    ISC correspond au sens singulet-triplet et RISC au sens triplet-singulet.
RADIATION : {NR : 5, fluorescent : 6}
    Dictionnaire contenant les modes de radiation. Les valeurs sont arbitraires.
BOUND_STATES : {bound : 7, unbound : 8}
    Dictionnaire contenant les états de liaison des excitons. Les valeurs sont arbitraires.

Fonctions
---------
is_electron(value : int) -> bool
    Vérifie si value correspond à la valeur de l'électron stockée dans le dictionnaire PARTICULES.
is_hole(value : int) -> bool
    Vérifie si value correspond à la valeur du trou stockée dans le dictionnaire PARTICULES.
is_exciton(value : int) -> bool
    Vérifie si value correspond à la valeur de l'exciton stockée dans le dictionnaire PARTICULES.
def is_ISC(value : int) -> bool
    Vérifie si value correspond à la valeur de ISC stockée dans le dictionnaire CROSSING.
is_RISC(value : int) -> bool
    Vérifie si value correspond à la valeur de RISC stockée dans le dictionnaire CROSSING.
is_NR(value : int) -> bool
    Vérifie si value correspond à la valeur de NR stockée dans le dictionnaire RADIATION.
is_fluorescent(value : int) -> bool
    Vérifie si value correspond à la valeur de fluorescent stockée dans le dictionnaire RADIATION.
is_bound(value : int) -> bool
    Vérifie si value correspond à la valeur de bound stockée dans le dictionnaire BOUND_STATES.
is_unbound(value : int) -> bool
    Vérifie si value correspond à la valeur de unbound stockée dans le dictionnaire BOUND_STATES.

Classes
-------
Point() : x, y, z
    Classe représentant un point. Les points peuvent s'additionner et se soutraire.
    Dans ce cas, le point est converti en vecteur, même dans le cas d'opération avec des nombres.
Vector(Point) : x, y, z
    Classe représentant un vecteur. Les vecteurs peuvent également se multiplier.
    La multiplication entre deux vecteurs donne le produit scalaire.
    La multiplication par un nombre donne un nouveau vecteur.

Event() : initial, final, tau
    Classe représentant un événement. Ceux-ci sont caractérisé par une position initiale,
    une position finale et une durée (tau).
Move(Event) : initial, final, tau, particule
    Classe représentant un événement de type déplacement. Contient aussi le type de particule à déplacer.
    Utilisé pour déplacer les charges au sein du réseau.
ISC(Event) : initial, final, tau, conversion
    Classe représentant un événement de type conversion intersystème. Contient aussi le sens de réaction.
    Utilisé pour les changement d'état de spin des excitons.
Decay(Event) : initial, final, tau, radiation
    Classe représentant un événement de type décroissance. Contient aussi le type de décroissance.
    Utilisé pour l'émission d'exciton.
Capture(Event) : initial, final, tau, particule
    Classe représentant un événement de type capture électronique. Contient aussi le type de particule.
    Utilisé pour la capture électronique aux électrodes opposées.
    Ne doit pas être utilisé pour des excitons.
Binding(Event), initial, final, tau, state
    Classe représentant un événement de type liaison/séparation. Contient aussi l'état de liaison.
    Utilisé pour la création et séparation des excitons. 
"""

from dataclasses import dataclass


EVENTS : dict[str, int] = {
    "move" : 1,
    "ISC" : 2,
    "decay" : 3,
    "capture" : 4,
    "unbound" : 5
}
PARTICULES : dict[str, int] = {
     "electron" : 1,
     "hole" : 2,
     "exciton" : 3
}


def is_electron(value : int) -> bool :
     return value == PARTICULES["electron"]

def is_hole(value : int) -> bool :
     return value == PARTICULES["hole"]

def is_exciton(value : int) -> bool :
     return value == PARTICULES["exciton"]


def is_move_event(value : int) -> bool :
    return value == EVENTS["move"]

def is_ISC_event(value : int) -> bool :
    return value == EVENTS["ISC"]

def is_decay_event(value : int) -> bool :
     return value == EVENTS["NR"]

def is_capture_event(value : int) -> bool :
     return value == EVENTS["capture"]

def is_unbound_event(value : int) -> bool :
     return value == EVENTS["unbound"]


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
        if issubclass(other, Point) :
            return Vector(self.x + other.x, self.y + other.y, self.z + other.z)
        elif isinstance(other, (float, int)) :
            return Vector(self.x + other, self.y + other, self.z + other)
        raise TypeError(f"other must be of type Vector, float or int, got {type(other)}")
            
    def __sub__(self, other) :
        if issubclass(other, Point) :
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
        if issubclass(other, Point) :
            return self.x * other.x + self.y * other.y + self.z * other.z
        elif isinstance(other, (float, int)) :
            return Vector(self.x * other, self.y * other, self.z * other)
        raise TypeError(f"other must be of type Vector, float or int, got {type(other)}")
    
    def __rmul__(self, other) :
        return self * other


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
            equality = [
                self.initial == other.initial,
                self.final == other.final,
                self.kind == other.kind,
                self.particule == other.particule
            ]
            return all(equality)
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