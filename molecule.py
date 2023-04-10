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

# Classe de type exciton utilisée par les molécules pour signaler la présence d'un exciton
class exciton :
    def __init__(self) :
        self.RNG = default_rng()
        self.spin = self._spin()    # spin de l'exciton

    # Choisi le spin de l'exciton : True pour singulet et False pour triplet
    def _spin(self) -> bool :
        random = self.RNG.uniform(0,100)
        if random < 25 :
            return True
        else :
            return False
# Classe de type energy utilisée pour stocker les niveaux d'énergie associé à chaque molécule
class energy :
    def __init__(self, HOMO : float, LUMO : float, S1: float, T1 : float) :
        self.HOMO = HOMO
        self.LUMO = LUMO
        self.S1 = S1
        self.T1 = T1
# Classe de type molecule, mère des classes fluorescent, tadf et host
class molecule :
    def __init__(self, position : point, voisins : list[point]) :
        self.RNG = default_rng()    # Seed de la molécule
        self.POS = position # Position de la molécule dans le réseau
        self.VOISINS = voisins
        self.electron = False   # Electron excité dans la LUMO
        self.hole = False   # Trou dans la HOMO
        self.exciton = None # Exciton (electron + trou)
        self.HOMO = None  # Energie de la HOMO (dépend du type de molécule)
        self.LUMO = None  # Energie de la LUMO (dépend du type de molécule)
        self.S1 = None  # Energie de l'état singulet (dépend du type de molécule)
        self.T1 = None  # Energie de l'état triplet (dépend du type de molécule)
    
    def _electron(self) :
        self.electron = not self.electron

    def _hole(self) :
        self.hole = not self.hole

    def _exciton(self) :
        if self.electron and self.hole :
            self.exciton = exciton()

    def _decay(self) :
        self.exciton = None
        self.electron = False
        self.hole = False
# Classe de type flourescent qui hérite des méthodes et des attribus de molecule. Utilisée dans le réseau.
class fluorescent(molecule) :
    def __init__(self, position : point, voisins : list[point], sigma : float = 0.1) :
        super().__init__(position, voisins)
        self.mu = energy(-5.3, -2.7, 2.69, 1.43)  # Moyenne pour le distribution gaussienne ! non implémenté
        self.sigma = sigma   # Déviation standard pour le distribution gaussienne ! non implémenté
        self._energie()

    def _energie(self) :
        self.HOMO = self.RNG.normal(self.mu.HOMO, self.sigma)
        self.LUMO = self.RNG.normal(self.mu.LUMO, self.sigma)
        self.S1 = self.RNG.normal(self.mu.S1, self.sigma)
        self.T1 = self.RNG.normal(self.mu.T1, self.sigma)
# Classe de type tadf qui hérite des méthodes et des attribus de molecule. Utilisée dans le réseau.
class tadf(molecule) :
    def __init__(self, position : point, voisins : list[point], sigma : float = 0.1) :
        super().__init__(position, voisins)
        self.mu = energy(-5.8, -2.6, 2.55, 2.52)  # Moyenne pour le distribution gaussienne ! non implémenté
        self.sigma = sigma   # Déviation standard pour le distribution gaussienne ! non implémenté
        self._energie()

    def _energie(self) :
        self.HOMO = self.RNG.normal(self.mu.HOMO, self.sigma)
        self.LUMO = self.RNG.normal(self.mu.LUMO, self.sigma)
        self.S1 = self.RNG.normal(self.mu.S1, self.sigma)
        self.T1 = self.RNG.normal(self.mu.T1, self.sigma)
# Classe de type host qui hérite des méthodes et des attribus de molecule. Utilisée dans le réseau.
class host(molecule) :
    def __init__(self, position : point, voisins : list[point], sigma : float = 0.1) :
        super().__init__(position, voisins)
        self.mu = energy(-6.0, -2.0, 3.50, 3.00)  # Moyenne pour le distribution gaussienne ! non implémenté
        self.sigma = sigma   # Déviation standard pour le distribution gaussienne ! non implémenté
        self._energie()

    def _energie(self) :
        self.HOMO = self.RNG.normal(self.mu.HOMO, self.sigma)
        self.LUMO = self.RNG.normal(self.mu.LUMO, self.sigma)
        self.S1 = self.RNG.normal(self.mu.S1, self.sigma)
        self.T1 = self.RNG.normal(self.mu.T1, self.sigma)
