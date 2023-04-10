#########################################################################################################
#
#   author(s) : Théo Piron
#   last update : 11/05/2022 (dd/mm/yyyy)
#   python version : 3.10.4
#   modules : reaseau, molecule, event, energy
#
#   Bugs connus :   - Fonction d'affichage pas efficace pour les grand réseaux
#
#   Remarques   :   Les énergies sont exprimées en (eV) et le temps en secondes
#
#########################################################################################################

# import matplotlib
from numpy import arange
# from TADF import tadf
# from emetteur import fluorescent
# from hote import host
from event import event, point
from molecule import exciton, fluorescent, tadf, host

from math import exp, log, floor
from numpy.random import default_rng
import matplotlib.pyplot as plt
from matplotlib import use


# Classe destinée à stocker les recombinaisons par type de molécule et dans l'ordre
class where :
    def __init__(self) :
        self.HOST = 0
        self.TADF = 0
        self.FLUO = 0
        self.Order = []
    
    def host(self) :
        self.HOST += 1
        self.Order.append(host)
    
    def tadf(self) :
        self.TADF += 1
        self.Order.append(tadf)

    def fluo(self) :
        self.FLUO += 1
        self.Order.append(fluorescent)

    def __str__(self) -> str:
        string = "Excitons formés sur les Host : " + str(self.HOST) + "\n"
        string = string + "Excitons formés sur les TADF : " + str(self.TADF) + "\n"
        string = string + "Excitons formés sur les Fluo : " + str(self.FLUO) + "\n"
        return string

    def __repr__(self) -> str:
        string = "Excitons formés sur les Host : " + str(self.HOST) + "\n"
        string = string + "Excitons formés sur les TADF : " + str(self.TADF) + "\n"
        string = string + "Excitons formés sur les Fluo : " + str(self.FLUO)
        return string

# Classe simulant un réseau cubique
class lattice :
    # Création d'un réseau tridimensionnel (sous forme d'une liste)
    # Arguments :   - dim : point contenant les dimensions du réseau
    #               - prop : tuple de taille 2 contenant les proportions (%) de matériaux tadf et fluo
    def __init__(self, dim : point, prop : tuple, Ez : float = 10**8, hashtag : int | None = None, sigma : float = 0.1) :
        # ID
        self.HASHTAG = hashtag
        # Vocabulaire du réseau
        self.REACTIONS = {"electron", "trou", "decay", "fluorescence"}   # Liste des événements possibles
        self.RNG = default_rng()    # Seed du réseau
        # Constantes du réseau
        self.DIM = dim  # dimensions du réseau
        self.PROP = prop    # pourcentages de matériaux TADF et Fluo
        self.E = point(0,0,Ez) # Champ électrique (0,0,Ez) (eV/m)
        self.A = 10**(-9)   # Constant de maille du réseau (m)
        self.NU = 10**13    # Fréquence de transfer de charge du réseau (Hz)
        self.KB = 8.617 * 10**(-5)  # Constante de Boltzmann (eV/K)
        self.T = 300    # Température du réseau (K)
        self.sigma = sigma  # Désordre énergétique (eV)
        self.COULOMB = 4.806*10**(-10) # Constante de Coulomb * charge élémentaire (Nm^2/C)
        # Réseau
        self.GRID : list[list[list[host | tadf | fluorescent]]] = [[[self._type(point(i,j,k)) for k in range(self.DIM.z)] for j in range(self.DIM.y)] for i in range(self.DIM.x)] # réseau
        # Evenements
        self.events_LUMO : list[event] = []   # stocke les événements liés à la LUMO (transfert de charge)
        self.events_HOMO: list[event] = []    # stocke les événements liés à la HOMO (transfert de charge)
        self.events_Exciton: list[event] = []  # stocke les événements d'émissions
        # Charges
        self.electrons = 10  # Nombre d'électrons
        self.electrons_pos : list[point] = []
        self.trous = 10  # Nombre de trous
        self.trous_pos : list[point] = []
        self.charges = 10
        self._injection()
        # IQE
        self.recombinaisons = 0 # Nombre de recombinaisons
        self.fluorescence = 0   # Nombre d'excitons fluorescents
        self.IQE = 0
        # Ordre
        self.where = where()    # Stocke le nombre de recombinaison par matériau
        # Temps
        self.temps = 0    # temps écoulé
    
    def __eq__(self, other) -> bool:
        if isinstance(other, lattice) :
            return self.HASHTAG == other.HASHTAG
        else :
            return False

    #################################################
    #
    #   Fonctions de démarrage du réseau
    #
    #################################################

    # Choisi le type de molécule de chaque noeud du réseau
    # Arguments :   - pos : position de la molécule dans le réseau
    def _type(self, pos : point) -> host | tadf | fluorescent :
        prop_tadf, prop_fluo = self.PROP
        random = self.RNG.uniform(0, 100)
        if random < 100 - prop_tadf - prop_fluo :
            return host(pos, self._calcul_voisins(pos), sigma = self.sigma)
        elif random < 100 - prop_fluo :
            return tadf(pos, self._calcul_voisins(pos), sigma = self.sigma)
        else :
            return fluorescent(pos, self._calcul_voisins(pos), sigma = self.sigma)

    # Détermine les voisins d'une molécule sur base d'un rayon d'action
    # Arguments :   - pos : position de la molécule dans le réseau
    def _calcul_voisins(self, pos : point) -> list[point] :
        voisins = []
        rayon = 1
        x_range = self._bvk_x(pos.x, rayon)
        y_range = self._bvk_y(pos.y, rayon)
        z_range = self._bvk_z(pos.z, rayon)
        for i in x_range :
            for j in y_range :
                for k in z_range :
                    voisins.append(point(i,j,k))
        voisins.remove(pos)
        return voisins

    # Détermine les indices x des voisins en tenant compte des conditions limites de Born von Kerman
    # Arguments :   - index : indice x
    #               - rayon : distance maximale
    def _bvk_x(self, index : int, rayon : int | float) -> list[int] :
        output = []
        for x in range(floor(index - rayon), floor(index + rayon) + 1) :
            if x < 0 :
                output.append(self.DIM.x + x)
                continue
            elif x > self.DIM.x - 1 :
                output.append(x - self.DIM.x)
                continue
            else :
                output.append(x)
                continue
        return output

    # Détermine les indices y des voisins en tenant compte des conditions limites de Born von Kerman
    # Arguments :   - index : indice y
    #               - rayon : distance maximale
    def _bvk_y(self, index : int, rayon : int | float) -> list[int] :
        output = []
        for y in range(floor(index - rayon), floor(index + rayon) + 1) :
            if y < 0 :
                output.append(self.DIM.y + y)
                continue
            elif y > self.DIM.y - 1 :
                output.append(y - self.DIM.y)
                continue
            else :
                output.append(y)
                continue
        return output

    # Détermine les indices z des voisins en tenant compte de la présence des électrodes.
    # Arguments :   - index : indice z
    #               - rayon : distance maximale
    def _bvk_z(self, index : int, rayon : int | float) -> list[int] :
        output = []
        for z in range(floor(index - rayon), floor(index + rayon) + 1) :
            if z < 0 :
                continue
            elif z > self.DIM.z - 1 :
                continue
            else :
                output.append(z)
                continue
        return output
    
    # Injecte les charges aux interfaces du réseau (z). Enregistre leurs positions et crée les premiers événements.
    def _injection(self) :
        while len(self.electrons_pos) < self.electrons :
            pos_x = self.RNG.integers(0, self.DIM.x)
            pos_y = self.RNG.integers(0, self.DIM.y)
            pos = point(pos_x, pos_y, self.DIM.z - 1)
            if pos not in self.electrons_pos :
                self.electrons_pos.append(pos)
            else :
                continue
        while len(self.trous_pos) < self.trous :
            pos_x = self.RNG.integers(0, self.DIM.x)
            pos_y = self.RNG.integers(0, self.DIM.y)
            pos = point(pos_x, pos_y, 0)
            if pos not in self.trous_pos :
                self.trous_pos.append(pos)
            else :
                continue
        for i in self.electrons_pos :
            self.GRID[i.x][i.y][i.z]._electron()
        for i in self.trous_pos :
            self.GRID[i.x][i.y][i.z]._hole()
        for i in self.electrons_pos + self.trous_pos :
            self._PPV(i)
        return

    ########################################################
    #
    #   Fonctions qui suppriment les événements obsolètes
    #
    ########################################################
    
    # Retire les événement obsolètes selon le type déterminé par First Reaction Method
    # Arguments :   - evenement : événement choisi par First Reaction Method sur lequel on travaille
    def _obsoletes(self, evenement : event) :
        reaction = evenement.REACTION
        # Désintégration d'un exciton ? 
        if reaction == "decay" or reaction == "fluorescence" :
            self._pop_decay(self.events_Exciton)
            return
        # Déplacement d'un électron ?
        elif reaction == "electron" :
            self._pop_electron(self.events_LUMO)
            # Si un exciton s'est formé, il faut nettoyer les trous également.
            if self._isExciton(evenement.final) :
                self._pop_hole(self.events_HOMO)
                return
            else :
                return
        # Déplacement d'un trou ?
        elif reaction == "trou" :
            self._pop_hole(self.events_HOMO)
            # Idem dans l'autre sens
            if self._isExciton(evenement.final) :
                self._pop_electron(self.events_LUMO)
                return
            else :
                return
    
    # Retire les événements de type exciton et decay
    # Arguments :   - tab : liste des événements sur lesquels on travaille (événements de type "exciton" et "decay")
    def _pop_decay(self, tab : list[event]) :
        obso = []
        for events in tab :
            # Molécules sans électrons
            if not self._isExciton(events.initial) :
                obso.append(events)
                continue
            continue
        for events in obso :
            tab.remove(events)
        return

    # Retire les événements de type déplacement d'électron
    # Arguments :   - tab : liste des événements sur lesquels on travaille (événements de type "electron")
    def _pop_electron(self, tab : list[event]) :
        obso = []
        for events in tab :
            # Molécules sans electrons
            if not self._isElectron(events.initial) :
                obso.append(events)
                continue
            # Molécules avec un exciton
            elif self._isExciton(events.initial) :
                obso.append(events)
                continue
            # Molécule avec un électron et dont le voisin est occupé
            elif self._isElectron(events.initial) and self._isElectron(events.final) :
                obso.append(events)
                continue
            continue
        for events in obso :
            tab.remove(events)
        return

    # Retire les événements de type déplacement de trou
    # Arguments :   - tab : liste des événements sur lesquels on travaille (événements de type "trou")
    def _pop_hole(self, tab : list[event]) :
        obso = []
        for events in tab :
            # Molécules vides
            if not self._isHole(events.initial) :
                obso.append(events)
                continue
            # Molécules avec un exciton
            elif self._isExciton(events.initial) :
                obso.append(events)
                continue
            # Molécule avec un trou et dont le voisin est occupé
            elif self._isHole(events.initial) and self._isHole(events.final) :
                obso.append(events)
                continue
            continue
        for events in obso :
            tab.remove(events)
        return
    
    ########################################################
    #
    #   Fonctions qui créent les nouveaux événements
    #
    ########################################################

    # Détermine le type de molécule sur lequel un exciton s'est formé
    # Arguments :   - pos : position à laquelle l'exciton s'est formé
    def _exciton_POS(self, pos : point) :
        if isinstance(self.GRID[pos.x][pos.y][pos.z], host) :
            self.where.host()
        elif isinstance(self.GRID[pos.x][pos.y][pos.z], tadf) :
            self.where.tadf()
        elif isinstance(self.GRID[pos.x][pos.y][pos.z], fluorescent) :
            self.where.fluo()
        else :
            print ("exciton_POS : exciton formé sur une molécule de type inconnue")

    # Vérifie le spin de l'exciton. Si singulet, ajoute l'événement
    # Arguments :   - pos : position de la molécule
    #               - react : type d'événement
    def _event_Exciton(self, pos : point) :
        self.recombinaisons += 1
        self._exciton_POS(pos)
        spin = self.GRID[pos.x][pos.y][pos.z].exciton.spin
        hote = isinstance(self.GRID[pos.x][pos.y][pos.z], host)
        # exciton non fluorescent (triplet ou sur hôte)
        if hote or not spin :
            react = "decay"
            evenement = event(pos, pos, 0., react)
            if evenement not in self.events_Exciton :
                self.events_Exciton.append(evenement)
            return
        # exciton fluorescent
        else :
            react = "fluorescence"
            evenement = event(pos, pos, 0., react)
            if evenement not in self.events_Exciton :
                self.events_Exciton.append(evenement)
            return

    # Calcul le temps de vie de l'exciton (inutilisée)
    def _Singulet(self, pos : point) -> float :
        dE = self.GRID[pos.x][pos.y][pos.z].S1
        return exp(dE)

    # Vérifie si la LUMO voisine est occupée. Sinon, ajoute l'événement si il n'est pas déjà présent dans la liste
    # Arguments :   - pos : position de la molécule
    #               - voisin : position cible
    #               - react : type d'événement
    def _event_LUMO(self, pos : point, voisin : point, react : str) :
        if not self._isElectron(voisin) :
            evenement = event(pos, voisin, self._LUMO(pos, voisin), react)
            if evenement not in self.events_LUMO :
                self.events_LUMO.append(evenement)
        return

    # Calcul le temps de réaction entre deux LUMO voisines.
    # Arguments :   - position, voisin : positions des molécules dont on calcule le temps de réaction
    def _LUMO(self, initial : point, final : point) -> float :
        rand = - self.RNG.uniform(-1,0)  # Nombre réel entre ]0,1]
        dist = self.A * (final - initial)
        dE = self.GRID[final.x][final.y][final.z].LUMO - self.GRID[initial.x][initial.y][initial.z].LUMO # LUMO
        dE += dist * self.E # Champ électrique
        dE += self._Coulomb_LUMO(initial, final) # Interaction Coulombienne
        if dE < 0 : 
            return - (log(rand)/self.NU)
        else :
            return - (log(rand)/self.NU) * exp(dE / (self.KB * self.T))
    
    # Calcul l'intéraction Coulombienne.
    # Arguments :   - initial : position de l'électron considéré
    def _Coulomb_LUMO(self, initial : point, final : point) -> float :
        output = 0
        for i in self.electrons_pos :
            if i == initial :
                continue
            elif not self._isExciton(i) :
                dist_init = self.A * (i - initial) ** 2
                dist_fin = self.A * (i - final) ** 2
                output += self.COULOMB * (1 / dist_fin - 1 / dist_init)
                continue
            else :
                continue
        for i in self.trous_pos :
            if i == final :
                output = -10**12
                return output
            elif not self._isExciton(i) :
                dist_init = self.A * (i - initial) ** 2
                dist_fin = self.A * (i - final) ** 2
                output -= self.COULOMB * (1 / dist_fin - 1 / dist_init)
                continue
            else :
                continue
        return output

    # Vérifie si la HOMO voisine est occupée. Sinon, ajoute l'événement si il n'est pas déjà présent dans la liste
    # Arguments :   - pos : position de la molécule
    #               - voisin : position cible
    #               - react : type d'événement
    def _event_HOMO(self, pos : point, voisin : point, react : str) :
        if not self._isHole(voisin) :
            evenement = event(pos, voisin, self._HOMO(pos, voisin), react)
            if evenement not in self.events_HOMO :
                self.events_HOMO.append(evenement)
        return

    # Calcul le temps de réaction entre deux HOMO voisines
    # Arguments :   - position, voisin : positions des molécules dont on calcule le temps de réaction
    def _HOMO(self, initial : point, final : point) -> float :
        rand = - self.RNG.uniform(-1,0)  # Nombre réel entre ]0,1]
        dist = self.A * (final - initial)
        dE = self.GRID[final.x][final.y][final.z].HOMO - self.GRID[initial.x][initial.y][initial.z].HOMO
        dE += dist * self.E
        dE = -dE
        dE += self._Coulomb_HOMO(initial, final)
        if dE < 0 : 
            return - (log(rand)/self.NU)
        else :
            return - (log(rand)/self.NU) * exp(dE / (self.KB * self.T))

    # Calcul l'intéraction Coulombienne.
    # Arguments :   - initial : position du trou considéré
    def _Coulomb_HOMO(self, initial : point, final : point) -> float :
        output = 0
        for i in self.trous_pos :
            if i == initial :
                continue
            elif not self._isExciton(i) :
                dist_init = self.A * (i - initial) ** 2
                dist_fin = self.A * (i - final) ** 2
                output += self.COULOMB * (1 / dist_fin - 1 / dist_init)
                continue
            else :
                continue
        for i in self.electrons_pos :
            if i == final :
                output = -10**12
                return output
            elif not self._isExciton(i) :
                dist_init = self.A * (i - initial) ** 2
                dist_fin = self.A * (i - final) ** 2
                output -= self.COULOMB * (1 / dist_fin - 1 / dist_init)
                continue
            else :
                continue
        return output

    # Détermine les événements d'une molécule selon le principe de plus proches voisins
    # Arguments :   - pos : position de la molécule sur laquelle on travaille
    def _PPV(self, pos : point) :
        # Si un exciton est présent sur la molécule, on crée un événement.
        if self._isExciton(pos) :
            self._event_Exciton(pos)
            return
        # Si un électron est présent sur la molécule, on fait le tour des voisins pour déterminer les événements.
        elif self._isElectron(pos) :
            react = "electron"
            for voisin in self.GRID[pos.x][pos.y][pos.z].VOISINS :
                self._event_LUMO(pos, voisin, react)
            return
        # Idem pour les trous
        elif self._isHole(pos) :
            react = "trou"
            for voisin in self.GRID[pos.x][pos.y][pos.z].VOISINS :
                self._event_HOMO(pos, voisin, react)
            return
        else :
            return
    
    # Détermine les événements des voisins d'une molécule anciennement occupée
    # Arguments :   - pos : position de la molécule sur laquelle on travaille
    def _old_PPV(self, pos : point) :
        for voisin in self.GRID[pos.x][pos.y][pos.z].VOISINS :
            # Si un exciton sur le voisin, on ne fait rien. (Géré par PPV)
            if self._isExciton(voisin) :
                continue
            # Si un électron sur le voisin, on crée un nouvel événement.
            elif self._isElectron(voisin) :
                react = "electron"
                self._event_LUMO(voisin, pos, react)
                continue
            # Idem pour un trou
            elif self._isHole(voisin) :
                react = "trou"
                self._event_HOMO(voisin, pos, react)
                continue
            else :
                continue
        return

    ########################################################
    #
    #   Fonctions qui déterminent la présence des charges
    #
    ########################################################

    # Vérifie la présence d'un électron
    # Arguments :   - point : position à laquelle on vérifie
    def _isElectron(self, point : point) -> bool :
        return self.GRID[point.x][point.y][point.z].electron

    # Vérifie la présence d'un trou
    # Arguments :   - point : position à laquelle on vérifie
    def _isHole(self, point : point) -> bool :
        return self.GRID[point.x][point.y][point.z].hole

    # Vérifie la présence d'un exciton 
    # Arguments :   - point : position à laquelle on vérifie   
    def _isExciton(self, point : point) -> bool :
        return isinstance(self.GRID[point.x][point.y][point.z].exciton, exciton)

    ####################################################
    #
    #   First Reaction Method
    #
    ####################################################

    # Détermine l'événement le plus rapide parmi ceux disponibles
    def _Fastest(self, tab : list[event]) -> event :
        tab.sort()
        return tab[0]

    # Algorithme de first reaction : localise l'événement qui prend le moins de temps, incrémente le temps, déplace la charge, vérifie la formation d'un exciton.
    def _First_Reaction_Method(self) -> bool :
        evenement : event | None = None
        # Cherche l'événement le plus rapide
        tab = self.events_HOMO + self.events_LUMO + self.events_Exciton
        if len(tab) > 0 :
            evenement = self._Fastest(tab)
        # Provoque l'événement
        if evenement == None :
            print ("First_Reaction : Aucun événement trouvé")
            return False
        elif evenement.REACTION in self.REACTIONS :
            # Déplace l'électron et vérifie la formation d'un excitons ou la capture par une électrode.
            if evenement.REACTION == "electron" :
                self.electrons_pos.remove(evenement.initial)
                self.GRID[evenement.initial.x][evenement.initial.y][evenement.initial.z]._electron()
                self.GRID[evenement.final.x][evenement.final.y][evenement.final.z]._electron()
                self.GRID[evenement.final.x][evenement.final.y][evenement.final.z]._exciton()
                if evenement.final.z == 0 and not self._isExciton(evenement.final) :
                    self.GRID[evenement.final.x][evenement.final.y][evenement.final.z]._electron()
                    self.electrons -= 1
                else :
                    self.electrons_pos.append(evenement.final)
            # Déplace le trou et vérifie la formation d'un excitons ou la capture par une électrode.
            elif evenement.REACTION == "trou" :
                self.trous_pos.remove(evenement.initial)
                self.GRID[evenement.initial.x][evenement.initial.y][evenement.initial.z]._hole()
                self.GRID[evenement.final.x][evenement.final.y][evenement.final.z]._hole()
                self.GRID[evenement.final.x][evenement.final.y][evenement.final.z]._exciton()
                if evenement.final.z == self.DIM.z-1 and not self._isExciton(evenement.final) :
                    self.GRID[evenement.final.x][evenement.final.y][evenement.final.z]._hole()
                    self.trous -= 1
                else :
                    self.trous_pos.append(evenement.final)
            # Désintègre l'exciton par fluorescence
            elif evenement.REACTION == "fluorescence" :
                self.fluorescence += 1
                self.electrons -= 1
                self.trous -= 1
                self.electrons_pos.remove(evenement.initial)
                self.trous_pos.remove(evenement.initial)
                self.GRID[evenement.initial.x][evenement.initial.y][evenement.initial.z]._decay()
            # Désintègre l'exciton de façon non radiative
            elif evenement.REACTION == "decay" :
                self.electrons -= 1
                self.trous -= 1
                self.electrons_pos.remove(evenement.initial)
                self.trous_pos.remove(evenement.initial)
                self.GRID[evenement.initial.x][evenement.initial.y][evenement.initial.z]._decay()
        else :
            print ("First_Reaction : Evénement de type inconnu")
            return False
        # Supprime les événements obsolètes
        self._obsoletes(evenement)
        # Fais passer le temps
        self.temps = self.temps + evenement.tau
        for i in self.events_HOMO + self.events_LUMO + self.events_Exciton :
            i.tau = i.tau - evenement.tau
        # Ajoute les nouvaux événements
        self._PPV(evenement.final)
        self._old_PPV(evenement.initial)
        return True

    # Fonction d'appel de l'algorithme de First Reaction Method. Permet de fixer un nombre d'itération maximal et s'arrête dans certains cas.
    # Arguments :   - stop :    Nombre d'itérations après lesquelles on arrête de faire évoluer le réseau
    #               - start :   Valeur de départ des itérations (par défaut 0)
    #                           Utile si on veut s'arrêter régulièrement pour afficher l'état du réseau et conserver le nombre d'étapes total.
    def Operations(self, stop : int, start : int = 0) -> None :
        for i in range(start, stop) :
            # Affiche l'étape régulièrement
            if not i%1000 :
                print ("Etape : " + str(i))
            # Plus de charge d'un type
            if self.electrons == 0 or self.trous == 0 :
                print("Opération : plus d'électrons ou de trous disponibles")
                print("Etape : " + str(i))
                if self.charges > 0 :
                    self.IQE = 100. * float(self.fluorescence)/float(self.charges)
                return
            temps = self.temps
            if self._First_Reaction_Method() :
                # Temps négatif : stop
                if temps > self.temps :
                    print ("Operations : Temps négatif détecté")
                    print ("Etape : " + str(i))
                    return
            # Pas d'événement trouvé : stop
            else :
                print ("Etape : " + str(i))
                print ("Fin")
                if self.charges > 0 :
                    self.IQE = 100. * float(self.fluorescence)/float(self.charges)
                return
        # S'assure que le réseau a eu le temps de recombiner des excitons
        if self.charges > 0 :
            self.IQE = 100. * float(self.fluorescence)/float(self.charges)
            return
        else :
            return
        
    ########################################################
    #
    #   Visualisation
    #
    ########################################################

    # Fonction d'affichage du réseau
    def Plot(self, nom : str) :
        use("Agg")
        if self.DIM.x*self.DIM.y*self.DIM.z > 50000 :
            print ("Plot : trop de molécules pour une représentation visuelle efficace")
            return
        plt.figure(dpi=100)
        axes = plt.axes(projection = "3d")
        # Création des points pour les molécules chargées
        electron : list[point] = []
        trou : list[point] = []
        exciton : list[point] = []
        legende = []
        for i in self.GRID :
            for j in i :
                for k in j :
                    if self._isExciton(k.POS) :
                        exciton.append(k.POS)
                        continue
                    elif self._isElectron(k.POS) :
                        electron.append(k.POS)
                        continue
                    elif self._isHole(k.POS) :
                        trou.append(k.POS)
                        continue
                    else :
                        continue
        # Création du nuage de point
        x, y, z = [], [], []
        if len(electron) > 0 :
            for i in electron :
                x.append(i.x)
                y.append(i.y)
                z.append(i.z)
            axes.scatter(x, y, z, marker = 'o', color = 'b')
            legende.append("electron")
            x.clear()
            y.clear()
            z.clear()
        if len(trou) > 0 :
            for i in trou :
                x.append(i.x)
                y.append(i.y)
                z.append(i.z)
            axes.scatter(x, y, z, marker = 'o', color = 'r')
            legende.append("hole")
            x.clear()
            y.clear()
            z.clear()
        if len(exciton) > 0 :
            for i in exciton :
                x.append(i.x)
                y.append(i.y)
                z.append(i.z)
            axes.scatter(x, y, z, marker = 'o', color = 'g')
            legende.append("exciton")
            x.clear()
            y.clear()
            z.clear()
        # Contour du réseau
        xGrid = [0, self.DIM.x-1]
        yGrid = [0, self.DIM.y-1]
        zGrid = [0, self.DIM.z-1]
        for i in xGrid :
            for j in yGrid :
                axes.plot([i,i], [j,j], [0, self.DIM.z-1], color = "k", linewidth = 1)
            for k in zGrid :
                axes.plot([i,i], [0, self.DIM.y-1], [k,k], color = "k", linewidth = 1)
        for j in yGrid :
            for k in zGrid :
                axes.plot([0, self.DIM.x-1], [j,j], [k,k], color = "k", linewidth = 1)
        # Options du graphique
        plt.legend(legende)
        axes.set_xlabel("x", size = 16)
        axes.set_ylabel("y", size = 16)
        axes.set_zlabel("z", size = 16)
        axes.set_xlim([0, self.DIM.x - 1])
        axes.set_ylim([0, self.DIM.y - 1])
        axes.set_zlim([0, self.DIM.z - 1])
        # Affichage
        plt.tight_layout()
        plt.savefig(nom, bbox_inches = "tight")
        plt.close()
        return

    ########################################################
    #
    #   Fonctions obsolètes
    #
    ########################################################

    # Détermine les voisins d'une molécule
    # Arguments :   - pos : position de la molécule dans le réseau
    def _voisins(self, pos : point) -> list[point] :
        voi = []
        voi.append(self._left(pos))
        voi.append(self._right(pos))
        voi.append(self._back(pos))
        voi.append(self._front(pos))
        if pos.z > 0 :
            voi.append(self._down(pos))
        if pos.z < self.DIM.z - 1 :
            voi.append(self._up(pos))
        return voi

    # Renvoie le voisin de gauche (x-1)
    # Argument :    - pos : position de la molécule
    def _left(self, pos : point) -> point :
        if pos.x == 0 :
            return point(self.DIM.x-1, pos.y, pos.z)
        else :
            return point(pos.x-1, pos.y, pos.z)
    
    # Renvoie le voisin de droite (x+1)
    # Argument :    - pos : position de la molécule
    def _right(self, pos : point) -> point :
        if pos.x == self.DIM.x - 1 :
            return point(0, pos.y, pos.z)
        else :
            return point(pos.x+1, pos.y, pos.z)

    # Renvoie le voisin de derrière (y-1)
    # Argument :    - pos : position de la molécule
    def _back(self, pos : point) -> point :
        if pos.y == 0 :
            return point(pos.x, self.DIM.y-1, pos.z)
        else :
            return point(pos.x, pos.y-1, pos.z)

    # Renvoie le voisin de devant (y+1)
    # Argument :    - pos : position de la molécule
    def _front(self, pos : point) -> point :
        if pos.y == self.DIM.y - 1 :
            return point(pos.x, 0, pos.z)
        else :
            return point(pos.x, pos.y+1, pos.z)

    # Renvoie le voisin du dessous (z-1)
    # Argument :    - pos : position de la molécule
    def _down(self, pos : point) -> point :
        return point(pos.x, pos.y, pos.z-1)

    # Renvoie le voisin du dessus (z+1)
    # Argument :    - pos : position de la molécule
    def _up(self, pos : point) -> point :
        return point(pos.x, pos.y, pos.z+1)

# test = lattice(point(5,5,5), (15,1))
# test.Operations(100)