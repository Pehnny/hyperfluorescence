#########################################################################################################
#
#   author(s) : Théo Piron
#   last update : 01/06/2023 (dd/mm/yyyy)
#   python version : 3.10.4
#   modules : reaseau, molecule, event, mu
#
#
#########################################################################################################
from numpy.random import default_rng, Generator
from event import point


class Exciton :
    """Classe représentant un exciton moléculaire.

    Attributes
    ----------
    spin : bool
        Etat de spin de l'exciton. True pour singulet et False pour triplet.
        
    Methods
    -------
    spin_state(random_number : int) -> bool
        Retourne l'état du spin de l'exciton lors de sa création.
    revert_spin() -> None
        Renverse la valeur de l'attribut spin.
    """

    def __init__(self, random_number : int) -> None :
        """Initialise l'instance de Exciton.

        Parameters
        ----------
        random_number : int
            Nombre aléatoire qui permet de déterminer l'état de spin de l'exciton.
        """
        self.spin : bool = self.spin_state(random_number)

    def __repr__(self) -> str:
        """Retourne l'état de spin de l'Exciton sous forme de texte.
        """
        return "singlet state exciton" if self.spin else "triplet state exciton"

    def spin_state(self, random_number : int) -> bool :
        """Génère l'état de spin de l'Exciton selon l'entrée.

        Parameters
        ----------
        random_number : int
            Nombre aléatoire qui permet de déterminer l'état de spin de l'exciton.

        Returns
        -------
        bool
            True si random_number est dans l'intervalle semi-ouvert [0,25). \n
            False si random_number est dans l'intervalle semi-ouvert [25,100).

        Raises
        ------
        ValueError
            Erreur levée quand le paramètre random_number est hors de l'intervalle
            semi-ouvert [0,100).
        """
        if 0 <= random_number < 25 :
            return True
        elif 25 <= random_number < 100 :
            return False
        raise ValueError(f"random_number = {random_number} but cannot exceed nor equal 100")
    
    def revert_spin(self) -> None :
        """Renverse l'état de spin de l'Exciton.
        """
        self.spin = not self.spin
        return


class Molecule :
    """Classe représentant une Molécule organique générique.

    Attributes
    ----------
    position : point
        Position de la molécule dans le réseau.
    neighbors : list[point]
        Liste des positions des voisins proches de la molécule.
    electron : bool
        Présence d'un électron dans la molécule.
    hole : bool
        Présence d'un trou dans la molécule.
    seed : Generator
        Graine de nombres pseudo-aléatoires propre à l'instance de la molécule.

    Methods
    ------- 
    switch_electron() -> None
        Renverse l'état de l'attribut electron.
    switch_electron() -> None
        Renverse l'état de l'attribut hole.
    generate_exciton(random_number : int) -> None
        Génère un exciton.
    exciton_decay() -> bool
        NotImplemented.
    """

    def __init__(self,
                 position : point,
                 voisins : list[point]) :
        """Initialise l'instance de Molecule.

        Parameters
        ----------
        position : point
            Position de la molécule dans le réseau.
        voisins : list[point]
            Liste des positions des voisins proches de la molécule.
        """
        self.position : point = position
        self.neighbors : list[point] = voisins
        self.electron : bool = False
        self.hole : bool = False
        self.seed : Generator = default_rng()
    
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
            self.exciton = Exciton(self.seed.integers(0,100))

    def exciton_decay(self) -> bool :
        """Détruit l'attribut exciton et remet les attributs electron et hole en False.
        """
        raise NotImplementedError(f"{__name__} is not implemented for {self.__class__.__name__}")


class Fluorescent(Molecule) :
    """Classe représentant une molécule fluorescente S1. Hérite de la classe Molecule.

    L'instance par défaut correspond à une molécume fluorescente bleue TBPe.

    Attributes
    ----------
    position : point
        Position de la molécule.
    neighbors : list[point]
        Liste des positions des voisins proches de la molécule.
    electron : bool
        Présence d'un électron dans la molécule.
    hole : bool
        Présence d'un trou dans la molécule.
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
    generate_electron(random_number : int) -> None
        Génère l'attribut exciton.
    exciton_decay() -> bool
        Décompose l'exciton. Remet les attributs electron et hole en False et supprime l'attribut
        exciton. Retourne True si l'exciton est singulet (émetteur fluorescent), False sinon.
    """

    def __init__(self, position : point,
                 voisins : list[point], homo_energy : float = -5.3,
                 lumo_energy : float = -2.7, s1_energy : float = 2.69,
                 t1_energy : float = 1.43, standard_deviation : float = 0.1) -> None :
        """Initialise l'instance de la classe Fluorescent.
        
        Parameters
        ----------
        position : point
            Position de la molécule dans le réseau.
        voisins : list[point]
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
        """Détruit l'attribut exciton et remet les attributs electron et hole en False.

        Returns
        -------
            Retourne True quand l'exciton est singulet (fluorescence), sinon False.

        Raises
        ------
        AttributeError
            Erreur levée quand la molécule n'a pas d'exciton moléculaire.
        """
        if hasattr(self, "exciton") :
            output : bool = self.exciton.spin
            self.electron = False
            self.hole = False
            del self.exciton
            return output
        else :
            raise AttributeError(f"Attribute exciton doesn't exist.")
        # try :
        #     output : bool = self.exciton.spin
        #     del self.exciton
        #     self.electron = False
        #     self.hole = False
        # except AttributeError :
        #     print("Warning ! User tried to decay exciton but no exciton was found.")
        #     output : bool = False
        # finally :
        #     return output


class TADF(Molecule) :
    """Classe représentant une molécule TADF S1. Hérite de la classe Molecule.

    L'instance par défaut correspond à une molécume fluorescente bleue ACRSA.

    Attributes
    ----------
    position : point
        Position de la molécule.
    neighbors : list[point]
        Liste des positions des voisins proches de la molécule.
    electron : bool
        Présence d'un électron dans la molécule.
    hole : bool
        Présence d'un trou dans la molécule.
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
    generate_electron(random_number : int) -> None
        Génère l'attribut exciton.
    exciton_decay() -> bool
        Décompose l'exciton. Remet les attributs electron et hole en False et supprime l'attribut
        exciton. Retourne True si l'exciton est singulet (émetteur fluorescent), False sinon.
    intersystem_crossing() -> None
        Converti le spin de l'exciton.
    """

    def __init__(self, position : point,
                 voisins : list[point], homo_energy : float = -5.8,
                 lumo_energy : float = -2.6, s1_energy : float = 2.55,
                 t1_energy : float = 2.52, standard_deviation : float = 0.1) -> None :
        """Initialise l'instance de la classe TADF.
        
        Parameters
        ----------
        position : point
            Position de la molécule dans le réseau.
        voisins : list[point]
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
        super().__init__(position, voisins)
        self.homo_energy : float = self.seed.normal(homo_energy, standard_deviation)
        self.lumo_energy : float = self.seed.normal(lumo_energy, standard_deviation)
        self.s1_energy : float = self.seed.normal(s1_energy, standard_deviation)
        self.t1_energy : float = self.seed.normal(t1_energy, standard_deviation)

    def exciton_decay(self) -> bool:
        """Détruit l'attribut exciton et remet les attributs electron et hole en False.

        Returns
        -------
            Retourne True quand l'exciton est singulet (fluorescence), sinon False.

        Raises
        ------
        AttributeError
            Erreur levée quand la molécule n'a pas d'exciton moléculaire.
        """
        if hasattr(self, "exciton") :
            output : bool = self.exciton.spin
            self.electron = False
            self.hole = False
            del self.exciton
            return output
        else :
            raise AttributeError(f"Attribute exciton doesn't exist.")
        
    def intersystem_crossing(self) -> None :
        """Converti l'état de spin de l'exciton.

        Raises
        ------
        AttributeError
            Erreur levée quand la molécule n'a pas d'exciton moléculaire.
        """
        if hasattr(self, "exciton") :
            self.exciton.revert_spin()
            return
        else :
            raise AttributeError("Attribute exciton doesn't exist.")


class Host(Molecule) :
    """Classe représentant une molécule Hote. Hérite de la classe Molecule.

    L'instance par défaut correspond à une molécume fluorescente bleue DPEPO.

    Attributes
    ----------
    position : point
        Position de la molécule.
    neighbors : list[point]
        Liste des positions des voisins proches de la molécule.
    electron : bool
        Présence d'un électron dans la molécule.
    hole : bool
        Présence d'un trou dans la molécule.
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
    generate_electron(random_number : int) -> None
        Génère l'attribut exciton.
    exciton_decay() -> bool
        Décompose l'exciton. Remet les attributs electron et hole en False et supprime l'attribut
        exciton. Retourne True si l'exciton est singulet (émetteur fluorescent), False sinon.
    """
    def __init__(self, position : point,
                 voisins : list[point], homo_energy : float = -6.0,
                 lumo_energy : float = -2.0, s1_energy : float = 3.50,
                 t1_energy : float = 3.00, standard_deviation : float = 0.1) -> None :
        """Initialise l'instance de la classe TADF.
        
        Parameters
        ----------
        position : point
            Position de la molécule dans le réseau.
        voisins : list[point]
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
        """Détruit l'attribut exciton et remet les attributs electron et hole en False.

        Returns
        -------
            Retourne False quel que soit l'état de spin car les molécules Hote n'émette pas dans
            le spectre du visible.

        Raises
        ------
        AttributeError
            Erreur levée quand la molécule n'a pas d'exciton moléculaire.
        """
        if hasattr(self, "exciton") :
            self.electron = False
            self.hole = False
            del self.exciton
            return False
        else :
            raise AttributeError(f"Attribute exciton doesn't exist.")