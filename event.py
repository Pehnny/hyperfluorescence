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


PARTICULES : dict[str, int] = {
    "electron" : 0,
    "hole" : 1,
    "exciton" : 2
}

CROSSING : dict[str, int] = {
    "ISC" : 3,
    "RISC" : 4
}

RADIATION : dict[str, int] = {
    "NR" : 5,
    "fluorescent" : 6
}

BOUND_STATE : dict[str, int] = {
    "bound" : 7,
    "unbound" : 8
}

def is_electron(value : int) -> bool :
        return value == PARTICULES["electron"]

def is_hole(value : int) -> bool :
        return value == PARTICULES["hole"]

def is_exciton(value : int) -> bool :
        return value == PARTICULES["exciton"]

def is_ISC(value : int) -> bool :
        return value == CROSSING["ISC"]

def is_RISC(value : int) -> bool :
        return value == CROSSING["RISC"]

def is_NR(value : int) -> bool :
        return value == RADIATION["NR"]

def is_fluorescent(value : int) -> bool :
        return value == RADIATION["fluorescent"]

def is_bound(value : int) -> bool :
        return value == BOUND_STATE["bound"]

def is_unbound(value : int) -> bool :
        return value == BOUND_STATE["unbound"]

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

    particule : int


@dataclass
class ISC(Event) :
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
        Sens de la conversion intersystème. Peut être "direct" (singlet-triplet) ou "reverse" (triplet-singlet). \n
        Les sens de réaction sont stockée dans la classe Crossing.
    """

    conversion : int


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

    radiation : int


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
        Type de particule sujette au déplacement. Les types autorisés sont "e-" ou "h+" \n
        Ces valeurs sont stockées comme attributs de classe dans la classe Particules.
    """

    particule : int


@dataclass
class Binding(Event) :
    """Classe représentant un événement de type séparation des charges d'un exciton.

    Opération de comparaison implémentée.

    Attributes
    ----------
    initial : Point
        Position de départ de l'événement.
    final : Point
        Position d'arrivée de l'événement.
    tau : float
        Durée de l'événement.
    bounding : str
        Type de particule sujette au déplacement. Les types autorisés sont "electron", "hole" et "exciton". \n
        Ces valeurs sont stockées comme attributs de classe dans la classe Particules.
    """

    state : int
