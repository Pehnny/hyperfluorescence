#########################################################################################################
#
#   author(s) : Théo Piron
#   last update : 01/06/2023 (dd/mm/yyyy)
#   python version : 3.10.4
#   modules : reaseau, molecule, event, mu
#
#   énergie des lumo définies négatives car occupé par des électrons virtuels (-e)
#   énergie des homo définies positives car occupé par des trous virtuels (+e)
#
#########################################################################################################
from event import Point
from random import Random
from dataclasses import dataclass
from math import inf

#   Exciton spin coupling state
EXCITON : dict[str, int] = {
    "none" : 0,
    "singlet" : 1,
    "doublet" : 2,
    "triplet" : 3
}

#   Transfer rates for different quantum mechanism [Hz].
#   NR : (triplet/singlet) non-radiative, F : fluorescence, PH : phosphorescence
TRANSFER_RATES : dict[str, float] = {
    "charges" : 10.**15,
    "DPEPO_NR" : inf,
    "ACRSA_F" : 4.58 * 10.**6,
    "ACRSA_PH" : 4.19 * 10.**6,
    "TBPe_F" : inf,
    "TBPe_NR" : inf
}

#   Förster energy transfer radius for ACRSA -> TBPe [nm].
#   STS : singlet to singlet, TTS : triplet to singlet.
TRANSFER_RADIUS : dict[str, float] = {
    "STS" : 5.55,
    "TTS" : 4.75
}

#   Spin Orbit Coupling amplitudes for ACRSA based on known values of ISC rates.
SOC : dict[str, float] = {
    "ISC" : 10.**8,
    "RISC" : 10.**5
}



class Molecule :
    """Classe représentant une Molécule organique générique.

    Attributes
    ----------
    position : Point
        Position de la molécule dans le réseau.
    neighbors : list[Point]
        Liste des positions des voisins proches de la molécule.
    electron : bool
        Présence d'un électron dans la molécule.
    hole : bool
        Présence d'un trou dans la molécule.
    has_exciton : bool
        Présence d'un exciton dans la molécule. Plus consistant avec electron et hole.
        Joue le rôle de hasattr().
    seed : Generator
        Graine de nombres pseudo-aléatoires propre à l'instance de la molécule.

    Methods
    ------- 
    switch_electron() -> None
        Renverse l'état de l'attribut electron.
    switch_electron() -> None
        Renverse l'état de l'attribut hole.
    generate_exciton() -> None
        Génère un exciton.
    unbound_exciton() -> None
        Sépare l'exciton en le remettant à 0.
    exciton_decay() -> bool
        NotImplemented.
    """

    def __init__(self, position : Point, neighbours : list[Point]) :
        """Initialise l'instance de Molecule.

        Parameters
        ----------
        position : Point
            Position de la molécule dans le réseau.
        voisins : list[Point]
            Liste des positions des voisins proches de la molécule.
        """
        self.position : Point = position
        self.neighbourhood : list[Point] = neighbours
        self.electron : bool = False
        self.hole : bool = False
        self.exciton : int = 0
        self.seed : Random = Random()
    
    def empty(self) -> bool :
        particules = [self.electron, self.hole, self.exciton]
        return not any(particules)
    
    def switch_electron(self) -> None :
        """Renverse l'état de l'attribut electron.
        """
        self.electron = not self.electron

    def switch_hole(self) -> None :
        """Renverse l'état de l'attribut hole.
        """
        self.hole = not self.hole

    def generate_exciton(self) -> None :       
        """Génère l'attribut exciton si les attributs electron et hole sont True.
        """
        if self.electron and self.hole :
            random = self.seed.random()
            if 0 <= random < 0.25 :
                self.exciton = EXCITON["singlet"]
            else :
                self.exciton = EXCITON["triplet"]

    def unbound_exciton(self) -> None :
        if self.exciton :
            self.exciton = EXCITON["none"]

    def exciton_decay(self) -> bool :
        """Méthode représentant la recombinaison d'un exciton, avec ou sans émission.
        """
        raise NotImplementedError(f"{self.exciton_decay.__name__} is not implemented.")



class Fluorophore(Molecule) :
    """Classe représentant une molécule fluorescente S1. Hérite de la classe Molecule.

    L'instance par défaut correspond à une molécume fluorescente bleue TBPe.

    Attributes
    ----------
    position : Point
        Position de la molécule.
    neighbors : list[Point]
        Liste des positions des voisins proches de la molécule.
    electron : bool
        Présence d'un électron dans la molécule.
    hole : bool
        Présence d'un trou dans la molécule.
    has_exciton : bool
        Présence d'un exciton dans la molécule. Plus consistant avec electron et hole.
        Joue le rôle de hasattr().
    seed : Generator
        Graine de nombres pseudo-aléatoires propre à l'instance de la molécule.
    homo_energy : float
        Energie de l'orbitale moléculaire occupée la plus haute générée aléatoirement selon
        une distrubition gaussienne dont la moyenne est l'argument éponyme.
    lumo_energy : float
        Energie de l'orbitale moléculaire inoccupée la plus basse générée aléatoirement selon
        une distrubition gaussienne dont la moyenne est l'argument éponyme.
    s1_energy : float
        Energie d'un exciton S1 au sein de la molécule générée aléatoirement selon
        une distrubition gaussienne dont la moyenne est l'argument éponyme.
    t1_energy : float
        Energie d'un exciton T1 au sein de la molécule générée aléatoirement selon
        une distrubition gaussienne dont la moyenne est l'argument éponyme.

    Methods
    -------
    switch_electron() -> None
        Renverse l'état de l'attribut electron.
    switch_electron() -> None
        Renverse l'état de l'attribut hole.
    generate_exciton() -> None
        Génère un exciton.
    unbound_exciton() -> None
        Sépare l'exciton en le remettant à 0.
    exciton_decay() -> bool
        Décompose l'exciton. Remet les attributs electron et hole en False et supprime l'attribut
        exciton. Retourne True si l'exciton est singulet (émetteur fluorescent), False sinon.
    """

    def __init__(self, position : Point, neighbours : list[Point],
                 homo_energy : float = 5.25, lumo_energy : float = -2.7, s1_energy : float = 2.69,
                 t1_energy : float = 1.43, standard_deviation : float = 0.1) -> None :
        """Initialise l'instance de la classe Fluorescent.
        
        Parameters
        ----------
        position : Point
            Position de la molécule dans le réseau.
        voisins : list[Point]
            Liste des positions des voisins proches de la molécule.
        homo_energy : float = -5.3 (5.25)
            Energie moyenne de l'orbitale moléculaire occupée la plus haute.
        lumo_energy : float = -2.7 (1.84)
            Energie moyenne de l'orbitale moléculaire inoccupée la plus basse.
        s1_energy : float = 2.69
            Energie moyenne du niveau d'énergie S1.
        t1_energy : float = 1.43
            Energie moyenne du niveau d'énergie T1.
        standard_deviation : float = 0.1
            Deviation standard des niveaux d'énergie.
        """
        super().__init__(position, neighbours)
        self.homo_energy : float = self.seed.gauss(homo_energy, standard_deviation)
        self.lumo_energy : float = self.seed.gauss(lumo_energy, standard_deviation)
        self.s1_energy : float = self.seed.gauss(s1_energy, standard_deviation)
        self.t1_energy : float = self.seed.gauss(t1_energy, standard_deviation)
        
    def exciton_decay(self) -> bool:
        """Méthode représentant la recombinaison d'un exciton, avec ou sans émission.

        Returns
        -------
            Si self.exciton != 0, change self.electron et self.hole à False, change self.exciton à 0.
            Enfin, retourne 1 si self.exciton était singulet car les Fluorescent sont fluorescentes, sinon 0.
        """

        if self.exciton :
            state = self.exciton
            self.exciton = EXCITON["none"]
            self.electron = False
            self.hole = False
            return state == EXCITON["singlet"]



class TADF(Molecule) :
    """Classe représentant une molécule TADF S1. Hérite de la classe Molecule.

    Les valeurs par défaut correspondent à la molécule TBPe.

    Attributes
    ----------
    position : Point
        Position de la molécule.
    neighbors : list[Point]
        Liste des positions des voisins proches de la molécule.
    electron : bool
        Présence d'un électron dans la molécule.
    hole : bool
        Présence d'un trou dans la molécule.
    has_exciton : bool
        Présence d'un exciton dans la molécule. Plus consistant avec electron et hole.
        Joue le rôle de hasattr().
    seed : Generator
        Graine de nombres pseudo-aléatoires propre à l'instance de la molécule.
    homo_energy : float
        Energie de l'orbitale moléculaire occupée la plus haute générée aléatoirement selon
        une distrubition gaussienne dont la moyenne est l'argument éponyme.
    lumo_energy : float
        Energie de l'orbitale moléculaire inoccupée la plus basse générée aléatoirement selon
        une distrubition gaussienne dont la moyenne est l'argument éponyme.
    s1_energy : float
        Energie d'un exciton S1 au sein de la molécule générée aléatoirement selon
        une distrubition gaussienne dont la moyenne est l'argument éponyme.
    t1_energy : float
        Energie d'un exciton T1 au sein de la molécule générée aléatoirement selon
        une distrubition gaussienne dont la moyenne est l'argument éponyme.

    Methods
    -------
    switch_electron() -> None
        Renverse l'état de l'attribut electron.
    switch_electron() -> None
        Renverse l'état de l'attribut hole.
    generate_exciton() -> None
        Génère un exciton.
    unbound_exciton() -> None
        Sépare l'exciton en le remettant à 0.
    exciton_decay() -> bool
        Décompose l'exciton. Remet les attributs electron et hole en False et supprime l'attribut
        exciton. Retourne True si l'exciton est singulet (émetteur fluorescent), False sinon.
    intersystem_crossing() -> None
        Converti le spin de l'exciton.
    """

    def __init__(self, position : Point,
                 neighbours : list[Point], homo_energy : float = 5.8,
                 lumo_energy : float = -2.6, s1_energy : float = 2.55,
                 t1_energy : float = 2.52, standard_deviation : float = 0.1) -> None :
        """Initialise l'instance de la classe TADF.

        Les valeurs part défaut correspondent à la molécule ACRSA
        
        Parameters
        ----------
        position : Point
            Position de la molécule dans le réseau.
        voisins : list[Point]
            Liste des positions des voisins proches de la molécule.
        homo_energy : float = -5.8
            Energie moyenne de l'orbitale moléculaire occupée la plus haute.
        lumo_energy : float = -2.6
            Energie moyenne de l'orbitale moléculaire inoccupée la plus basse.
        s1_energy : float = 2.55
            Energie moyenne du niveau d'énergie S1.
        t1_energy : float = 2.52
            Energie moyenne du niveau d'énergie T1.
        standard_deviation : float = 0.1
            Deviation standard des niveaux d'énergie.
        """
        super().__init__(position, neighbours)
        self.homo_energy : float = self.seed.gauss(homo_energy, standard_deviation)
        self.lumo_energy : float = self.seed.gauss(lumo_energy, standard_deviation)
        self.s1_energy : float = self.seed.gauss(s1_energy, standard_deviation)
        self.t1_energy : float = self.seed.gauss(t1_energy, standard_deviation)

    def exciton_decay(self) -> bool:
        """Méthode représentant la recombinaison d'un exciton, avec ou sans émission.

        Returns
        -------
            Si self.exciton != 0, change self.electron et self.hole à False, change self.exciton à 0.
            Enfin, retourne 1 si self.exciton était singulet car les TADF sont fluorescentes, sinon 0.
        """

        if self.exciton :
            self.exciton = EXCITON["none"]
            self.electron = False
            self.hole = False
            return False
        
    def intersystem_crossing(self) -> None :
        """Converti l'état de spin de l'exciton.

        Si self.exciton est un sigulet, self.exciton est changé en triplet et vice versa.
        """
        if self.exciton == EXCITON["singlet"] :
            self.exciton = EXCITON["triplet"]
        elif self.exciton == EXCITON["triplet"] :
            self.exciton = EXCITON["singlet"]



class Host(Molecule) :
    """Classe représentant une molécule Hote. Hérite de la classe Molecule.

    Les valeurs par défaut correspondent à la molécule DPEPO.

    Attributes
    ----------
    position : Point
        Position de la molécule.
    neighbors : list[Point]
        Liste des positions des voisins proches de la molécule.
    electron : bool
        Présence d'un électron dans la molécule.
    hole : bool
        Présence d'un trou dans la molécule.
    has_exciton : bool
        Présence d'un exciton dans la molécule. Plus consistant avec electron et hole.
        Joue le rôle de hasattr().
    seed : Generator
        Graine de nombres pseudo-aléatoires propre à l'instance de la molécule.
    homo_energy : float
        Energie de l'orbitale moléculaire occupée la plus haute générée aléatoirement selon
        une distrubition gaussienne dont la moyenne est l'argument éponyme.
    lumo_energy : float
        Energie de l'orbitale moléculaire inoccupée la plus basse générée aléatoirement selon
        une distrubition gaussienne dont la moyenne est l'argument éponyme.
    s1_energy : float
        Energie d'un exciton S1 au sein de la molécule générée aléatoirement selon
        une distrubition gaussienne dont la moyenne est l'argument éponyme.
    t1_energy : float
        Energie d'un exciton T1 au sein de la molécule générée aléatoirement selon
        une distrubition gaussienne dont la moyenne est l'argument éponyme.

    Methods
    -------
    switch_electron() -> None
        Renverse l'état de l'attribut electron.
    switch_electron() -> None
        Renverse l'état de l'attribut hole.
    generate_exciton() -> None
        Génère un exciton.
    unbound_exciton() -> None
        Sépare l'exciton en le remettant à 0.
    exciton_decay() -> bool
        Décompose l'exciton. Remet les attributs electron et hole en False et supprime l'attribut
        exciton. Retourne True si l'exciton est singulet (émetteur fluorescent), False sinon.
    """

    def __init__(self, position : Point, neighbours : list[Point],
                 homo_energy : float = 6.0, lumo_energy : float = -2.0, s1_energy : float = 3.50,
                 t1_energy : float = 3.00, standard_deviation : float = 0.1) -> None :
        """Initialise l'instance de la classe Host.
        
        Parameters
        ----------
        position : Point
            Position de la molécule dans le réseau.
        voisins : list[Point]
            Liste des positions des voisins proches de la molécule.
        homo_energy : float = -6.0
            Energie moyenne de l'orbitale moléculaire occupée la plus haute.
        lumo_energy : float = -2.0
            Energie moyenne de l'orbitale moléculaire inoccupée la plus basse.
        s1_energy : float = 3.50
            Energie moyenne du niveau d'énergie S1.
        t1_energy : float = 3.00
            Energie moyenne du niveau d'énergie T1.
        standard_deviation : float = 0.1
            Deviation standard des niveaux d'énergie.
        """
        super().__init__(position, neighbours)
        self.homo_energy : float = self.seed.gauss(homo_energy, standard_deviation)
        self.lumo_energy : float = self.seed.gauss(lumo_energy, standard_deviation)
        self.s1_energy : float = self.seed.gauss(s1_energy, standard_deviation)
        self.t1_energy : float = self.seed.gauss(t1_energy, standard_deviation)
        
    def exciton_decay(self) -> bool:
        """Méthode représentant la recombinaison d'un exciton, avec ou sans émission.

        Returns
        -------
            Si self.exciton != 0, change self.electron et self.hole à False, change self.exciton à 0.
            Enfin, retourne 0 car les Host ne sont pas des molécules émittrice (dans le visible).
        """

        if self.exciton :
            self.exciton = EXCITON["none"]
            self.electron = False
            self.hole = False
            return False


@dataclass
class Proportion :
    """Classe représentant les proportions de chaque molécules au sein du réseau

    Attirbutes
    ----------
    host : float
        Proportion de molécules Host au sein du réseau.
    tadf : float
        Proportion de molécules TADF au sein du réseau.
    fluo : float
        Proportion de molécules Fluorescent au sein du réseau.
    """
    host : float
    tadf : float
    fluo : float
