#########################################################################################################
#
#   author(s) : Théo Piron
#   last update : 11/05/2022 (dd/mm/yyyy)
#   python version : 3.10.4
#   modules : reaseau, molecule, event, mu
#
#
#########################################################################################################

from numpy.random import default_rng, Generator
from event import point


class Exciton :
    """
    Class Exciton : (random_number : int)
        Classe de type exciton moléculaire.

        Arguments :
            random_number : int
                Nombre aléatoire qui permet de déterminer l'état de spin de l'exciton.
        Instance attributs :
            spin : bool
                Etat de spin de l'exciton. True pour singulet et False pour triplet
        Methodes : 
            __spin : (random_number : int) -> bool
                Retourne l'état du spin de l'exciton.
    """  
    def __init__(self,
                 random_number : int) -> None :
        self.spin : bool = self.__spin(random_number)

    def __repr__(self) -> str:
        if self.spin :
            return "singlet state exciton"
        else :
            return "triplet state exciton"

    def __spin(self,
               random_number : int) -> bool :
        if random_number < 25 :
            return True
        elif random_number < 100 :
            return False
        else :
            raise ValueError(f"random_number = {random_number} but cannot exceed nor equal 100")


class Molecule :
    """
    Class Molecule : (position : point, voisins : list[point])
        Classe de base pour toutes les molécules.

        Arguments :
            position : point
                Position de la molécule.
            voisins : list[point]
                Liste des positions des voisins proches de la molécule.
        Instance attributs :
            position : point
                Position de la molécule.
            neighbors : list[point]
                Liste des positions des voisins proches de la molécule.
            electron : bool
                Présence d'un électron dans la molécule.
            hole : bool
                Présence d'un trou dans la molécule.
            exciton : Exciton
                Exciton moléculaire de type singulet ou triplet. Attribut temporaire.
        Methodes : 
            switch_electron : () -> None
                Renverse l'état de l'attribut electron.
            switch_electron : () -> None
                Renverse l'état de l'attribut hole.
            generate_electron : (random_number : int) -> None
                Génère un exciton.
            exciton_decay : () -> bool
                NotImplemented.
    """  
    def __init__(self,
                 position : point,
                 voisins : list[point]) :
        self.position : point = position
        self.neighbors : list[point] = voisins
        self.electron : bool = False
        self.hole : bool = False
        self.seed : Generator = default_rng()
    
    def switch_electron(self) -> None :
        self.electron = not self.electron

    def switch_hole(self) -> None :
        self.hole = not self.hole

    def generate_exciton(self) -> None :
        if self.electron and self.hole :
            self.exciton = Exciton(self.seed.integers())

    def exciton_decay(self) -> bool :     
        raise NotImplementedError(f"{__name__} is not implemented for {self.__class__.__name__}")


class Fluorescent(Molecule):
    """
    Class Fluorescent : (position : point, voisins : list[point], homo_energy : float = -5.3,
                        lumo_energy : float = -2.7, s1_energy : float = 2.69, t1_energy : float = 1.43,
                        standard_deviation : float = 0.1)
        Sous-classe de Molecule. 
        Les valeurs par défaut correspondent à la molécule fluorescente TBPe.

        Super arguments :
            position : point
                Position de la molécule.
            voisins : list[point]
                Liste des positions des voisins proches de la molécule.
        Arguments :
            homo_energy : float = -5.3
                Energie moyenne de l'orbitale moléculaire occupée la plus haute.
            lumo_energy : float = -2.7
                Energie moyenne de l'orbitale moléculaire inoccupée la plus basse.
            s1_energy : float = 2.69
                Energie moyenne d'un exciton S1 au sein de la molécule.
            t1_energy : float = 1.43
                Energie moyenne d'un exciton T1 au sein de la molécule.
            standard_deviation : float = 0.1
                Déviation standard des niveaux d'énergie de la molécule.
        Super instance attributs :
            position : point
                Position de la molécule.
            neighbors : list[point]
                Liste des positions des voisins proches de la molécule.
            electron : bool
                Présence d'un électron dans la molécule.
            hole : bool
                Présence d'un trou dans la molécule.
            exciton : Exciton
                Exciton moléculaire de type singulet ou triplet. Attribut temporaire.
        Instance attributs :
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
        Super methods :
            switch_electron : () -> None
                Renverse l'état de l'attribut electron.
            switch_electron : () -> None
                Renverse l'état de l'attribut hole.
            generate_electron : (random_number : int) -> None
                Génère un exciton.
        Methods :
            exciton_decay : () -> bool
                Décompose l'exciton. Remet les attributs electron et hole en False et supprime l'attribut
                exciton. Retourne True si l'exciton est singulet (émetteur fluorescent), False sinon.
    """
    def __init__(self,
                 position : point,
                 voisins : list[point],
                 homo_energy : float = -5.3,
                 lumo_energy : float = -2.7,
                 s1_energy : float = 2.69,
                 t1_energy : float = 1.43,
                 standard_deviation : float = 0.1) -> None :
        super().__init__(position,
                         voisins)
        self.homo_energy : float = self.seed.normal(homo_energy,
                                                    standard_deviation)
        self.lumo_energy : float = self.seed.normal(lumo_energy,
                                                    standard_deviation)
        self.s1_energy : float = self.seed.normal(s1_energy,
                                                  standard_deviation)
        self.t1_energy : float = self.seed.normal(t1_energy,
                                                  standard_deviation)
        
    def exciton_decay(self) -> bool:
        try :
            output : bool = self.exciton.spin
            del self.exciton
            self.electron = False
            self.hole = False
        except AttributeError :
            print("Warning ! User tried to decay exciton but no exciton was found.")
            output : bool = False
        finally :
            return output


class TADF(Molecule) :
    """
    Class TADF :
        Sous-classe de Molecule. 
        Les valeurs par défaut correspondent à la molécule TADF de type ACRSA.

        Super arguments :
            position : point
                Position de la molécule.
            voisins : list[point]
                Liste des positions des voisins proches de la molécule.
        Arguments :
            homo_energy : float = -5.3
                Energie moyenne de l'orbitale moléculaire occupée la plus haute.
            lumo_energy : float = -2.7
                Energie moyenne de l'orbitale moléculaire inoccupée la plus basse.
            s1_energy : float = 2.69
                Energie moyenne d'un exciton S1 au sein de la molécule.
            t1_energy : float = 1.43
                Energie moyenne d'un exciton T1 au sein de la molécule.
            standard_deviation : float = 0.1
                Déviation standard des niveaux d'énergie de la molécule.
        Super instance attributs :
            position : point
                Position de la molécule.
            neighbors : list[point]
                Liste des positions des voisins proches de la molécule.
            electron : bool
                Présence d'un électron dans la molécule.
            hole : bool
                Présence d'un trou dans la molécule.
            exciton : Exciton
                Exciton moléculaire de type singulet ou triplet. Attribut temporaire.
        Instance attributs :
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
        Super methods :
            switch_electron : () -> None
                Renverse l'état de l'attribut electron.
            switch_electron : () -> None
                Renverse l'état de l'attribut hole.
            generate_electron : (random_number : int) -> None
                Génère un exciton.
        Methods :
            exciton_decay : () -> bool
                Décompose l'exciton. Remet les attributs electron et hole en False et supprime l'attribut
                exciton. Retourne True si l'exciton est singulet (émetteur fluorescent), False sinon.
            intersystem_crossing : () -> None
                Renverse l'état de spin de l'exciton.
    """
    def __init__(self,
                 position : point,
                 voisins : list[point],
                 homo_energy : float = -5.8,
                 lumo_energy : float = -2.6,
                 s1_energy : float = 2.55,
                 t1_energy : float = 2.52,
                 standard_deviation : float = 0.1) -> None :
        super().__init__(position,
                         voisins)
        self.homo_energy : float = self.seed.normal(homo_energy,
                                                    standard_deviation)
        self.lumo_energy : float = self.seed.normal(lumo_energy,
                                                    standard_deviation)
        self.s1_energy : float = self.seed.normal(s1_energy,
                                                  standard_deviation)
        self.t1_energy : float = self.seed.normal(t1_energy,
                                                  standard_deviation)

    def exciton_decay(self) -> bool:
        try :
            output : bool = self.exciton.spin
            del self.exciton
            self.electron = False
            self.hole = False
        except AttributeError :
            print("Warning ! User tried to decay exciton but no exciton was found.")
            output : bool = False
        finally :
            return output
        
    def intersystem_crossing(self) -> None :
        try :
            self.exciton.spin = not self.exciton.spin
        except AttributeError :
            print("Warning ! User tried to change exciton spin state but no exciton was found.")
        finally :
            return

class Host(Molecule) :
    """
    Class Host :
        Sous-classe de Molecule. 
        Les valeurs par défaut correspondent à la molécule hôte de type DPEPO.

        Super arguments :
            position : point
                Position de la molécule.
            voisins : list[point]
                Liste des positions des voisins proches de la molécule.
        Arguments :
            homo_energy : float = -5.3
                Energie moyenne de l'orbitale moléculaire occupée la plus haute.
            lumo_energy : float = -2.7
                Energie moyenne de l'orbitale moléculaire inoccupée la plus basse.
            s1_energy : float = 2.69
                Energie moyenne d'un exciton S1 au sein de la molécule.
            t1_energy : float = 1.43
                Energie moyenne d'un exciton T1 au sein de la molécule.
            standard_deviation : float = 0.1
                Déviation standard des niveaux d'énergie de la molécule.
        Super instance attributs :
            position : point
                Position de la molécule.
            neighbors : list[point]
                Liste des positions des voisins proches de la molécule.
            electron : bool
                Présence d'un électron dans la molécule.
            hole : bool
                Présence d'un trou dans la molécule.
            exciton : Exciton
                Exciton moléculaire de type singulet ou triplet. Attribut temporaire.
        Instance attributs :
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
        Super methods :
            switch_electron : () -> None
                Renverse l'état de l'attribut electron.
            switch_electron : () -> None
                Renverse l'état de l'attribut hole.
            generate_electron : (random_number : int) -> None
                Génère un exciton.
        Methods :
            exciton_decay : () -> bool
                Décompose l'exciton. Remet les attributs electron et hole en False et supprime l'attribut
                exciton. Retourne False (non émetteur).
    """
    def __init__(self,
                 position : point,
                 voisins : list[point],
                 homo_energy : float = -6.0,
                 lumo_energy : float = -2.0,
                 s1_energy : float = 3.50,
                 t1_energy : float = 3.00,
                 standard_deviation : float = 0.1) -> None :
        super().__init__(position,
                         voisins)
        self.homo_energy : float = self.seed.normal(homo_energy,
                                                    standard_deviation)
        self.lumo_energy : float = self.seed.normal(lumo_energy,
                                                    standard_deviation)
        self.s1_energy : float = self.seed.normal(s1_energy,
                                                  standard_deviation)
        self.t1_energy : float = self.seed.normal(t1_energy,
                                                  standard_deviation)
        
    def exciton_decay(self) -> bool:
        try :
            del self.exciton
            self.electron = False
            self.hole = False
        except AttributeError :
            print("Warning ! User tried to decay exciton but no exciton was found.")
        finally :
            return False