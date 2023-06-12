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
from event import *
from molecule import *
from constants import *

from math import exp, log, floor
import matplotlib.pyplot as plt
from matplotlib import use
from random import Random

# Classe destinée à stocker les recombinaisons par type de molécule et dans l'ordre
# class where :
#     def __init__(self) :
#         self.HOST = 0
#         self.TADF = 0
#         self.FLUO = 0
#         self.Order = []
    
#     def host(self) :
#         self.HOST += 1
#         self.Order.append(host)
    
#     def tadf(self) :
#         self.TADF += 1
#         self.Order.append(tadf)

#     def fluo(self) :
#         self.FLUO += 1
#         self.Order.append(fluorescent)

#     def __str__(self) -> str:
#         string = "Excitons formés sur les Host : " + str(self.HOST) + "\n"
#         string = string + "Excitons formés sur les TADF : " + str(self.TADF) + "\n"
#         string = string + "Excitons formés sur les Fluo : " + str(self.FLUO) + "\n"
#         return string

#     def __repr__(self) -> str:
#         string = "Excitons formés sur les Host : " + str(self.HOST) + "\n"
#         string = string + "Excitons formés sur les TADF : " + str(self.TADF) + "\n"
#         string = string + "Excitons formés sur les Fluo : " + str(self.FLUO)
#         return string


class lattice :
    def __init__(self, dimension : Point,
                 proportion : Proportion, Ez : float = 10.**8,
                 charges : int = 10) -> None :
        # Vocabulaire du réseau
        self.SEED : Random = Random()
        # Constantes du réseau
        self.DIMENSION : Point =  dimension
        self.PROPORTION : Proportion = proportion
        self.ELECTRIC_FIELD : Point = Point(0,0,Ez)     # [eV/m]
        self.LATTICE_CONSTANT : float = 1.              # [nm]
        self.TRANSFER_RATE : float = 10.**13            # [Hz]
        self.TEMPERATURE : float = 300.                 # [K]
        # Réseau
        self.GRID : list[list[list[Host | TADF | Fluorescent]]] = self._construction()
        # Evenements
        self.events_LUMO : list[Event] = []   # stocke les événements liés à la LUMO (transfert de charge)
        self.events_HOMO: list[Event] = []    # stocke les événements liés à la HOMO (transfert de charge)
        self.events_Exciton: list[Event] = []  # stocke les événements d'émissions
        # Charges
        self.charges : int = charges  # Nombre de charges de chaque signe (!= charges totales)
        self.positions_electrons : list[Point] = self._injection(self.DIMENSION.z - 1, charges)
        self.positions_holes : list[Point] = self._injection(0, charges)
        for electron, hole in zip(self.positions_electrons, self.positions_holes):
            self.GRID[electron.x][electron.y][electron.z].switch_electron()
            # self._PPV(electron)
            self.GRID[hole.x][hole.y][hole.z].switch_hole()
            # self._PPV(trou)
        # IQE
        self.recombinaisons = 0 # Nombre de recombinaisons
        self.fluorescence = 0   # Nombre d'excitons fluorescents
        self.iqe = 0    # Efficacité quantique interne
        # Ordre
        self.where = where()    # Stocke le nombre de recombinaison par matériau
        # Temps
        self.temps : float = 0.    # temps écoulé

    def _construction(self) -> list[list[list[Host | TADF | Fluorescent]]] :
        if self.DIMENSION.z < 3 :
            raise ValueError(f"The z component of dimension must be > 2, got {self.DIMENSION.z}.")
        x_max : int = self.DIMENSION.x
        y_max : int = self.DIMENSION.y
        z_max : int = self.DIMENSION.z
        grid_size : int = x_max * y_max * z_max
        n_fluo : int = int(grid_size * self.PROPORTION.fluo)
        n_tadf : int = int(grid_size * self.PROPORTION.tadf)
        n_host : int = grid_size - 2 * x_max * y_max - n_fluo - n_tadf
        sub_z_max : int = z_max - 2
        sub_grid_size : int = x_max * y_max * sub_z_max

        sub_grid : list[int] = self.SEED.sample(
            [0, 1, 2],
            k = n_host + n_tadf + n_fluo,
            counts = [n_host, n_tadf, n_fluo]
        )

        if len(sub_grid) != sub_grid_size :
            raise IndexError(f"Size of sub_grid ({len(sub_grid)}) and z_max*y_max*x_max ({sub_grid_size}) must match")
        
        boundaries : list[list[int]] = [[0 for j in range(y_max)] for i in range(x_max)]
        grid : list[list[list[int]]] = [
            boundaries,
            [[sub_grid[i + x_max * j + x_max * y_max * k] for i in range(x_max)] for j in range(y_max)] for k in range(sub_z_max)
        ]
        grid.append(boundaries)
        
        for x in range(x_max):
            for y in range(y_max):
                for z in range(z_max):
                    value = grid[x,y,z]
                    position = Point(x,y,z)
                    grid[x,y,z] == self._molecule_type(value, position)
        return grid
    
    def _molecule_type(self, n : int, pos : Point) :
        if n == 0 :
            return Host(pos, self._neighbourhood(pos))
        elif n == 1 :
            return TADF(pos, self._neighbourhood(pos))
        elif n == 2 :
            return Fluorescent(pos, self._neighbourhood(pos))
        raise ValueError(f"n should be 0, 1 or 2, got {n}")

    def _neighbourhood(self, pos : Point, distance : int = 1) -> list[Point] :
        x_range = self._born_von_karman(pos.x, distance, "x")
        y_range = self._born_von_karman(pos.y, distance, "y")
        z_range = self._born_von_karman(pos.z, distance, "z")
        return [Point(x,y,z) for x in x_range for y in y_range for z in z_range if (x, y, z) != (pos.x, pos.y, pos.z)]

    def _born_von_karman(self, position : int, distance : int, axe : str) -> list[int] :
        size = getattr(self.DIMENSION, axe)
        lower_bound = position - distance
        upper_bound = position + distance
        if lower_bound > -1 and upper_bound < size :
            return list(range(position - distance, position + distance + 1))
        elif lower_bound < 0 :
            if axe != "z" :
                return list(range(0, distance + 1)) + list(range(size - distance, size))
            else :
                return list(range(0, distance + 1))
        else :
            if axe != "z" :
                return list(range(position - distance, size)) + list(range(distance))
            else :
                return list(range(position - distance, size))
    
    def _injection(self, z : int, charges : int) -> list[Point] :
        molecules = self.DIMENSION.x * self.DIMENSION.y
        if charges > molecules :
            raise ValueError(f"Required {charges} charges but only {molecules} molecules available.")
        possibilities = (Point(x, y, z) for x in range(self.DIMENSION.x) for y in self.DIMENSION.y)
        return self.SEED.sample(possibilities, k = charges)

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
        rand = - self.rng.uniform(-1,0)  # Nombre réel entre ]0,1]
        dist = self.LATTICE_CONSTANT * (final - initial)
        dE = self.GRID[final.x][final.y][final.z].LUMO - self.GRID[initial.x][initial.y][initial.z].LUMO # LUMO
        dE += dist * self.ELECTRIC_FIELD # Champ électrique
        dE += self._Coulomb_LUMO(initial, final) # Interaction Coulombienne
        if dE < 0 : 
            return - (log(rand)/self.TRANSFER_RATE)
        else :
            return - (log(rand)/self.TRANSFER_RATE) * exp(dE / (self.constanteBoltzmann * self.TEMPERATURE))
    
    # Calcul l'intéraction Coulombienne.
    # Arguments :   - initial : position de l'électron considéré
    def _Coulomb_LUMO(self, initial : point, final : point) -> float :
        output = 0
        for i in self.electronPositions :
            if i == initial :
                continue
            elif not self._isExciton(i) :
                dist_init = self.LATTICE_CONSTANT * (i - initial) ** 2
                dist_fin = self.LATTICE_CONSTANT * (i - final) ** 2
                output += self.coulomb * (1 / dist_fin - 1 / dist_init)
                continue
            else :
                continue
        for i in self.trouPositions :
            if i == final :
                output = -10**12
                return output
            elif not self._isExciton(i) :
                dist_init = self.LATTICE_CONSTANT * (i - initial) ** 2
                dist_fin = self.LATTICE_CONSTANT * (i - final) ** 2
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
        dist = self.LATTICE_CONSTANT * (final - initial)
        dE = self.GRID[final.x][final.y][final.z].HOMO - self.GRID[initial.x][initial.y][initial.z].HOMO
        dE += dist * self.ELECTRIC_FIELD
        dE = -dE
        dE += self._Coulomb_HOMO(initial, final)
        if dE < 0 : 
            return - (log(rand)/self.TRANSFER_RATE)
        else :
            return - (log(rand)/self.TRANSFER_RATE) * exp(dE / (self.constanteBoltzmann * self.TEMPERATURE))

    # Calcul l'intéraction Coulombienne.
    # Arguments :   - initial : position du trou considéré
    def _Coulomb_HOMO(self, initial : point, final : point) -> float :
        output = 0
        for i in self.trouPositions :
            if i == initial :
                continue
            elif not self._isExciton(i) :
                dist_init = self.LATTICE_CONSTANT * (i - initial) ** 2
                dist_fin = self.LATTICE_CONSTANT * (i - final) ** 2
                output += self.coulomb * (1 / dist_fin - 1 / dist_init)
                continue
            else :
                continue
        for i in self.electronPositions :
            if i == final :
                output = -10**12
                return output
            elif not self._isExciton(i) :
                dist_init = self.LATTICE_CONSTANT * (i - initial) ** 2
                dist_fin = self.LATTICE_CONSTANT * (i - final) ** 2
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
        elif evenement.REACTION in self.reactions :
            # Déplace l'électron et vérifie la formation d'un excitons ou la capture par une électrode.
            if evenement.REACTION == "electron" :
                self.electronPositions.remove(evenement.initial)
                self.GRID[evenement.initial.x][evenement.initial.y][evenement.initial.z].Electron()
                self.GRID[evenement.final.x][evenement.final.y][evenement.final.z].Electron()
                self.GRID[evenement.final.x][evenement.final.y][evenement.final.z]._exciton()
                if evenement.final.z == 0 and not self._isExciton(evenement.final) :
                    self.GRID[evenement.final.x][evenement.final.y][evenement.final.z].Electron()
                    self.electrons -= 1
                else :
                    self.electronPositions.append(evenement.final)
            # Déplace le trou et vérifie la formation d'un excitons ou la capture par une électrode.
            elif evenement.REACTION == "trou" :
                self.trouPositions.remove(evenement.initial)
                self.GRID[evenement.initial.x][evenement.initial.y][evenement.initial.z].Trou()
                self.GRID[evenement.final.x][evenement.final.y][evenement.final.z].Trou()
                self.GRID[evenement.final.x][evenement.final.y][evenement.final.z]._exciton()
                if evenement.final.z == self.DIMENSION.z-1 and not self._isExciton(evenement.final) :
                    self.GRID[evenement.final.x][evenement.final.y][evenement.final.z].Trou()
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
                self.GRID[evenement.initial.x][evenement.initial.y][evenement.initial.z]._decay()
            # Désintègre l'exciton de façon non radiative
            elif evenement.REACTION == "decay" :
                self.electrons -= 1
                self.trous -= 1
                self.electronPositions.remove(evenement.initial)
                self.trouPositions.remove(evenement.initial)
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
        xgrid = [0, self.DIMENSION[0]-1]
        ygrid = [0, self.DIMENSION[1]-1]
        zgrid = [0, self.DIMENSION[2]-1]
        for i in xgrid :
            for j in ygrid :
                axes.plot([i,i], [j,j], [0, self.DIMENSION[2]-1], color = "k", linewidth = 1)
            for k in zgrid :
                axes.plot([i,i], [0, self.DIMENSION[1]-1], [k,k], color = "k", linewidth = 1)
        for j in ygrid :
            for k in zgrid :
                axes.plot([0, self.DIMENSION[0]-1], [j,j], [k,k], color = "k", linewidth = 1)
        # Options du graphique
        plt.legend(legende)
        axes.set_xlabel("x", size = 16)
        axes.set_ylabel("y", size = 16)
        axes.set_zlabel("z", size = 16)
        axes.set_xlim([0, self.DIMENSION[0] - 1])
        axes.set_ylim([0, self.DIMENSION[1] - 1])
        axes.set_zlim([0, self.DIMENSION[2] - 1])
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
        prop_tadf, prop_fluo = self.PROPORTION
        random = self.rng.uniform(0, 100)
        if random < 100 - prop_tadf - prop_fluo :
            return host(pos, self.__neighbourhood(pos), sigma = self.sigma)
        elif random < 100 - prop_fluo :
            return tadf(pos, self.__neighbourhood(pos), sigma = self.sigma)
        else :
            return fluorescent(pos, self.__neighbourhood(pos), sigma = self.sigma)

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
        if pos.z < self.DIMENSION.z - 1 :
            voi.append(self._up(pos))
        return voi

    # Renvoie le voisin de gauche (x-1)
    # Argument :    - pos : position de la molécule
    def _left(self, pos : point) -> point :
        if pos.x == 0 :
            return point(self.DIMENSION.x-1, pos.y, pos.z)
        else :
            return point(pos.x-1, pos.y, pos.z)
    
    # Renvoie le voisin de droite (x+1)
    # Argument :    - pos : position de la molécule
    def _right(self, pos : point) -> point :
        if pos.x == self.DIMENSION.x - 1 :
            return point(0, pos.y, pos.z)
        else :
            return point(pos.x+1, pos.y, pos.z)

    # Renvoie le voisin de derrière (y-1)
    # Argument :    - pos : position de la molécule
    def _back(self, pos : point) -> point :
        if pos.y == 0 :
            return point(pos.x, self.DIMENSION.y-1, pos.z)
        else :
            return point(pos.x, pos.y-1, pos.z)

    # Renvoie le voisin de devant (y+1)
    # Argument :    - pos : position de la molécule
    def _front(self, pos : point) -> point :
        if pos.y == self.DIMENSION.y - 1 :
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