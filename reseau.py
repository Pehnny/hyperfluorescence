#########################################################################################################
#
#   author(s) : Théo Piron
#   last update : 01/02/2023 (dd/mm/yyyy)
#   python version : 3.10.4
#   numpy version : 1.21.4
#   matplotlib version : 3.5.0
#   modules : reseau, molecule, event, energy
#
#   Bugs connus :   - Fonction d'affichage pas efficace pour les grand réseaux
#
#   Remarques   :   Les énergies sont exprimées en (eV) et le temps en secondes
#
#########################################################################################################

import numpy as np
from event import event, point
from molecule import exciton, fluorescent, tadf, host

from math import exp, log, floor
from numpy.random import default_rng
import matplotlib.pyplot as plt
from matplotlib import use

from tools import elapsedTime

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
    @elapsedTime("Init :")
    def __init__(self, dimension : tuple, proportions : tuple, Ez : float = 10**8, hashtag : int | None = None, sigma : float = 0.1) :
        # ID
        self.HASHTAG = hashtag
        # Vocabulaire du réseau
        self.reactions = {"electron", "trou", "decay", "fluorescence"}   # Liste des événements possibles
        self.rng = default_rng(12345)    # Seed du réseau
        # Constantes du réseau
        self.dimension =  dimension # dimensions du réseau
        self.proportion = proportions    # pourcentages de matériaux TADF et Fluo
        self.electriField = point(0,0,Ez) # Champ électrique (0,0,Ez) (eV/m)
        self.constanteMaille = 10**(-9)   # Constant de maille du réseau (m)
        self.frequenceTransfert = 10**13    # Fréquence de transfer de charge du réseau (Hz)
        self.constanteBoltzmann = 8.617 * 10**(-5)  # Constante de Boltzmann (eV/K)
        self.temperature = 300    # Température du réseau (K)
        self.sigma = sigma  # Désordre énergétique (eV)
        self.coulomb = 4.806*10**(-10) # Constante de Coulomb * charge élémentaire (Nm^2/C)
        # Réseau
        self.grid = self.__Construction()
        # Evenements
        self.events_LUMO : list[event] = []   # stocke les événements liés à la LUMO (transfert de charge)
        self.events_HOMO: list[event] = []    # stocke les événements liés à la HOMO (transfert de charge)
        self.events_Exciton: list[event] = []  # stocke les événements d'émissions
        # Charges
        self.charges = 10  # Nombre de charges de chaque signe (!= charges totales)
        self.electronPositions = self.__Injection(self.dimension[2]-1)
        self.trouPositions = self.__Injection(0)
        for electron, trou in zip(self.electronPositions, self.trouPositions):
            self.grid[electron.x][electron.y][electron.z].Electron()
            # self._PPV(electron)
            self.grid[trou.x][trou.y][trou.z].Trou()
            # self._PPV(trou)
        # IQE
        self.recombinaisons = 0 # Nombre de recombinaisons
        self.fluorescence = 0   # Nombre d'excitons fluorescents
        self.iqe = 0    # Efficacité quantique interne
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
    def __Construction(self) -> list[list[list[host | tadf | fluorescent]]] :
        taille = self.dimension[0]*self.dimension[1]*self.dimension[2]
        propTADF, propFluo = self.proportion
        propTADF, propFluo = floor(taille * propTADF / 100), floor(taille * propFluo / 100)
        valeur = np.zeros(taille, int)
        for i in range(propTADF):
            valeur[i] = 1
        for i in range(propTADF, propTADF+propFluo):
            valeur[i] = 2
        self.rng.shuffle(valeur)
        valeur = valeur.reshape((self.dimension[0], self.dimension[1], self.dimension[2]))
        valeur = valeur.tolist()
        for x in range(self.dimension[0]):
            for y in range(self.dimension[1]):
                for z in range(self.dimension[2]):
                    pos = point(x,y,z)
                    valeur[x][y][z] = self.__MoleculeType(valeur[x][y][z], pos)
        return valeur
    
    def __MoleculeType(self, n : int, pos : point) :
        if n == 0 :
            return host(pos, self.__CalculVoisins(pos), sigma=self.sigma)
        elif n == 1 :
            return tadf(pos, self.__CalculVoisins(pos), sigma=self.sigma)
        else :
            return fluorescent(pos, self.__CalculVoisins(pos), sigma=self.sigma)

    # Détermine les voisins d'une molécule sur base d'un rayon d'action
    # Arguments :   - pos : position de la molécule dans le réseau
    def __CalculVoisins(self, pos : point) -> list[point] :
        distance = 1
        xRange = self.__BVKX(pos.x, distance)
        yRange = self.__BVKY(pos.y, distance)
        zRange = self.__BVKZ(pos.z, distance)
        return [point(x, y, z) for x in xRange for y in yRange for z in zRange if (x, y, z) != (pos.x, pos.y, pos.z)]

    # Détermine les indices x des voisins en tenant compte des conditions limites de Born von Kerman
    # Arguments :   - pos : indice x à la positions évaluée
    #               - distance : demi arrête du cube autour de la position x
    def __BVKX(self, pos : int, distance : int) -> list[int] :
        borneInf = pos - distance
        borneSup = pos + distance
        if borneInf > -1 and borneSup < self.dimension[0] :
            return list(range(pos-distance, pos+distance+1))
        elif borneInf < 0 :
            return list(range(0, distance+1)) + list(range(self.dimension[0]-distance, self.dimension[0]))
        else :
            return list(range(pos-distance, self.dimension[0])) + list(range(distance))

    # Détermine les indices y des voisins en tenant compte des conditions limites de Born von Kerman
    # Arguments :   - pos : indice y à la positions évaluée
    #               - distance : demi arrête du cube autour de la position y
    def __BVKY(self, pos : int, distance : int) -> list[int] :
        borneInf = pos - distance
        borneSup = pos + distance
        if borneInf > -1 and borneSup < self.dimension[1] :
            return list(range(pos-distance, pos+distance+1))
        elif borneInf < 0 :
            return list(range(0, distance+1)) + list(range(self.dimension[1]-distance, self.dimension[1]))
        else :
            return list(range(pos-distance, self.dimension[1])) + list(range(distance))

    # Détermine les indices z des voisins en tenant compte de la présence des électrodes.
    # Arguments :   - pos : indice z à la positions évaluée
    #               - distance : demi arrête du cube autour de la position z
    def __BVKZ(self, pos : int, distance : int | float) -> list[int] :
        borneInf = pos - distance
        borneSup = pos + distance
        if borneInf > -1 and borneSup < self.dimension[2] :
            return list(range(pos-distance, pos+distance+1))
        elif borneInf < 0 :
            return list(range(0, distance+1))
        else :
            return list(range(pos-distance, self.dimension[2]))
    
    # Injecte les charges aux interfaces du réseau (z). Enregistre leurs positions et crée les premiers événements.
    @elapsedTime("Injection :")
    def __Injection(self, z : int) -> list[point] :
        posX = self.rng.choice(self.dimension[0], size = self.charges, replace=False, shuffle=False)
        posY = self.rng.choice(self.dimension[1], size = self.charges, replace=False, shuffle=False)
        return [point(x, y, z) for x,y in zip(posX, posY)]

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
        if isinstance(self.grid[pos.x][pos.y][pos.z], host) :
            self.where.host()
        elif isinstance(self.grid[pos.x][pos.y][pos.z], tadf) :
            self.where.tadf()
        elif isinstance(self.grid[pos.x][pos.y][pos.z], fluorescent) :
            self.where.fluo()
        else :
            print ("exciton_POS : exciton formé sur une molécule de type inconnue")

    # Vérifie le spin de l'exciton. Si singulet, ajoute l'événement
    # Arguments :   - pos : position de la molécule
    #               - react : type d'événement
    def _event_Exciton(self, pos : point) :
        self.recombinaisons += 1
        self._exciton_POS(pos)
        spin = self.grid[pos.x][pos.y][pos.z].exciton.spin
        hote = isinstance(self.grid[pos.x][pos.y][pos.z], host)
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
        dE = self.grid[pos.x][pos.y][pos.z].S1
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
        rand = - self.rng.uniform(-1,0)  # Nombre réel entre ]0,1]
        dist = self.constanteMaille * (final - initial)
        dE = self.grid[final.x][final.y][final.z].LUMO - self.grid[initial.x][initial.y][initial.z].LUMO # LUMO
        dE += dist * self.electriField # Champ électrique
        dE += self._Coulomb_LUMO(initial, final) # Interaction Coulombienne
        if dE < 0 : 
            return - (log(rand)/self.frequenceTransfert)
        else :
            return - (log(rand)/self.frequenceTransfert) * exp(dE / (self.constanteBoltzmann * self.temperature))
    
    # Calcul l'intéraction Coulombienne.
    # Arguments :   - initial : position de l'électron considéré
    def _Coulomb_LUMO(self, initial : point, final : point) -> float :
        output = 0
        for i in self.electronPositions :
            if i == initial :
                continue
            elif not self._isExciton(i) :
                dist_init = self.constanteMaille * (i - initial) ** 2
                dist_fin = self.constanteMaille * (i - final) ** 2
                output += self.coulomb * (1 / dist_fin - 1 / dist_init)
                continue
            else :
                continue
        for i in self.trouPositions :
            if i == final :
                output = -10**12
                return output
            elif not self._isExciton(i) :
                dist_init = self.constanteMaille * (i - initial) ** 2
                dist_fin = self.constanteMaille * (i - final) ** 2
                output -= self.coulomb * (1 / dist_fin - 1 / dist_init)
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
        rand = - self.rng.uniform(-1,0)  # Nombre réel entre ]0,1]
        dist = self.constanteMaille * (final - initial)
        dE = self.grid[final.x][final.y][final.z].HOMO - self.grid[initial.x][initial.y][initial.z].HOMO
        dE += dist * self.electriField
        dE = -dE
        dE += self._Coulomb_HOMO(initial, final)
        if dE < 0 : 
            return - (log(rand)/self.frequenceTransfert)
        else :
            return - (log(rand)/self.frequenceTransfert) * exp(dE / (self.constanteBoltzmann * self.temperature))

    # Calcul l'intéraction Coulombienne.
    # Arguments :   - initial : position du trou considéré
    def _Coulomb_HOMO(self, initial : point, final : point) -> float :
        output = 0
        for i in self.trouPositions :
            if i == initial :
                continue
            elif not self._isExciton(i) :
                dist_init = self.constanteMaille * (i - initial) ** 2
                dist_fin = self.constanteMaille * (i - final) ** 2
                output += self.coulomb * (1 / dist_fin - 1 / dist_init)
                continue
            else :
                continue
        for i in self.electronPositions :
            if i == final :
                output = -10**12
                return output
            elif not self._isExciton(i) :
                dist_init = self.constanteMaille * (i - initial) ** 2
                dist_fin = self.constanteMaille * (i - final) ** 2
                output -= self.coulomb * (1 / dist_fin - 1 / dist_init)
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
            for voisin in self.grid[pos.x][pos.y][pos.z].VOISINS :
                self._event_LUMO(pos, voisin, react)
            return
        # Idem pour les trous
        elif self._isHole(pos) :
            react = "trou"
            for voisin in self.grid[pos.x][pos.y][pos.z].VOISINS :
                self._event_HOMO(pos, voisin, react)
            return
        else :
            return
    
    # Détermine les événements des voisins d'une molécule anciennement occupée
    # Arguments :   - pos : position de la molécule sur laquelle on travaille
    def _old_PPV(self, pos : point) :
        for voisin in self.grid[pos.x][pos.y][pos.z].VOISINS :
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
        return self.grid[point.x][point.y][point.z].electron

    # Vérifie la présence d'un trou
    # Arguments :   - point : position à laquelle on vérifie
    def _isHole(self, point : point) -> bool :
        return self.grid[point.x][point.y][point.z].hole

    # Vérifie la présence d'un exciton 
    # Arguments :   - point : position à laquelle on vérifie   
    def _isExciton(self, point : point) -> bool :
        return isinstance(self.grid[point.x][point.y][point.z].exciton, exciton)

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
        elif evenement.REACTION in self.reactions :
            # Déplace l'électron et vérifie la formation d'un excitons ou la capture par une électrode.
            if evenement.REACTION == "electron" :
                self.electronPositions.remove(evenement.initial)
                self.grid[evenement.initial.x][evenement.initial.y][evenement.initial.z].Electron()
                self.grid[evenement.final.x][evenement.final.y][evenement.final.z].Electron()
                self.grid[evenement.final.x][evenement.final.y][evenement.final.z]._exciton()
                if evenement.final.z == 0 and not self._isExciton(evenement.final) :
                    self.grid[evenement.final.x][evenement.final.y][evenement.final.z].Electron()
                    self.electrons -= 1
                else :
                    self.electronPositions.append(evenement.final)
            # Déplace le trou et vérifie la formation d'un excitons ou la capture par une électrode.
            elif evenement.REACTION == "trou" :
                self.trouPositions.remove(evenement.initial)
                self.grid[evenement.initial.x][evenement.initial.y][evenement.initial.z].Trou()
                self.grid[evenement.final.x][evenement.final.y][evenement.final.z].Trou()
                self.grid[evenement.final.x][evenement.final.y][evenement.final.z]._exciton()
                if evenement.final.z == self.dimension.z-1 and not self._isExciton(evenement.final) :
                    self.grid[evenement.final.x][evenement.final.y][evenement.final.z].Trou()
                    self.trous -= 1
                else :
                    self.trouPositions.append(evenement.final)
            # Désintègre l'exciton par fluorescence
            elif evenement.REACTION == "fluorescence" :
                self.fluorescence += 1
                self.electrons -= 1
                self.trous -= 1
                self.electronPositions.remove(evenement.initial)
                self.trouPositions.remove(evenement.initial)
                self.grid[evenement.initial.x][evenement.initial.y][evenement.initial.z]._decay()
            # Désintègre l'exciton de façon non radiative
            elif evenement.REACTION == "decay" :
                self.electrons -= 1
                self.trous -= 1
                self.electronPositions.remove(evenement.initial)
                self.trouPositions.remove(evenement.initial)
                self.grid[evenement.initial.x][evenement.initial.y][evenement.initial.z]._decay()
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
                    self.iqe = 100. * float(self.fluorescence)/float(self.charges)
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
                    self.iqe = 100. * float(self.fluorescence)/float(self.charges)
                return
        # S'assure que le réseau a eu le temps de recombiner des excitons
        if self.charges > 0 :
            self.iqe = 100. * float(self.fluorescence)/float(self.charges)
            return
        else :
            return
        
    ########################################################
    #
    #   Visualisation
    #
    ########################################################
    # Fonction d'affichage du réseau
    @elapsedTime("Plot :")
    def Plot(self, nom : str = "Inconnu") :
        use("Agg")
        plt.figure(dpi=100)
        axes = plt.axes(projection = "3d")
        # Création des points pour les molécules chargées
        electrons = self.electronPositions
        trous = self.trouPositions
        excitons = [electron for electron, trou in zip(electrons, trous) if electron == trou]
        legende = []
        # Création du nuage de point
        x, y, z = [], [], []
        if len(electrons) > 0 :
            for pos in electrons :
                x.append(pos.x)
                y.append(pos.y)
                z.append(pos.z)
            axes.scatter(x, y, z, marker = 'o', color = 'b')
            legende.append("electrons")
            x.clear()
            y.clear()
            z.clear()
        if len(trous) > 0 :
            for pos in trous :
                x.append(pos.x)
                y.append(pos.y)
                z.append(pos.z)
            axes.scatter(x, y, z, marker = 'o', color = 'r')
            legende.append("trous")
            x.clear()
            y.clear()
            z.clear()
        if len(excitons) > 0 :
            for pos in excitons :
                x.append(pos.x)
                y.append(pos.y)
                z.append(pos.z)
            axes.scatter(x, y, z, marker = 'o', color = 'g')
            legende.append("exciton")
            x.clear()
            y.clear()
            z.clear()
        # Contour du réseau
        xgrid = [0, self.dimension[0]-1]
        ygrid = [0, self.dimension[1]-1]
        zgrid = [0, self.dimension[2]-1]
        for i in xgrid :
            for j in ygrid :
                axes.plot([i,i], [j,j], [0, self.dimension[2]-1], color = "k", linewidth = 1)
            for k in zgrid :
                axes.plot([i,i], [0, self.dimension[1]-1], [k,k], color = "k", linewidth = 1)
        for j in ygrid :
            for k in zgrid :
                axes.plot([0, self.dimension[0]-1], [j,j], [k,k], color = "k", linewidth = 1)
        # Options du graphique
        plt.legend(legende)
        axes.set_xlabel("x", size = 16)
        axes.set_ylabel("y", size = 16)
        axes.set_zlabel("z", size = 16)
        axes.set_xlim([0, self.dimension[0] - 1])
        axes.set_ylim([0, self.dimension[1] - 1])
        axes.set_zlim([0, self.dimension[2] - 1])
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
    # Choisi le type de molécule de chaque noeud du réseau
    # Arguments :   - pos : position de la molécule dans le réseau
    def __TypeOld(self, pos : point) -> host | tadf | fluorescent :
        prop_tadf, prop_fluo = self.proportion
        random = self.rng.uniform(0, 100)
        if random < 100 - prop_tadf - prop_fluo :
            return host(pos, self.__CalculVoisins(pos), sigma = self.sigma)
        elif random < 100 - prop_fluo :
            return tadf(pos, self.__CalculVoisins(pos), sigma = self.sigma)
        else :
            return fluorescent(pos, self.__CalculVoisins(pos), sigma = self.sigma)

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
        if pos.z < self.dimension.z - 1 :
            voi.append(self._up(pos))
        return voi

    # Renvoie le voisin de gauche (x-1)
    # Argument :    - pos : position de la molécule
    def _left(self, pos : point) -> point :
        if pos.x == 0 :
            return point(self.dimension.x-1, pos.y, pos.z)
        else :
            return point(pos.x-1, pos.y, pos.z)
    
    # Renvoie le voisin de droite (x+1)
    # Argument :    - pos : position de la molécule
    def _right(self, pos : point) -> point :
        if pos.x == self.dimension.x - 1 :
            return point(0, pos.y, pos.z)
        else :
            return point(pos.x+1, pos.y, pos.z)

    # Renvoie le voisin de derrière (y-1)
    # Argument :    - pos : position de la molécule
    def _back(self, pos : point) -> point :
        if pos.y == 0 :
            return point(pos.x, self.dimension.y-1, pos.z)
        else :
            return point(pos.x, pos.y-1, pos.z)

    # Renvoie le voisin de devant (y+1)
    # Argument :    - pos : position de la molécule
    def _front(self, pos : point) -> point :
        if pos.y == self.dimension.y - 1 :
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

test = lattice((50,50,50), (15,1))
test.Plot("test")
# test.Operations(100)