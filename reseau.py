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
import constants as cst

from math import exp, log
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
    """Classe représentant un réseau cristallin de type OLED hyperfluorescente.

    lattice(dimension : tuple[int,int,int], proportion : tuple[float,float,float],
            electric_field : float = 10.**8, charges : int = 10)

    Attributes
    ----------
    _seed : Random
        Graine de nombres pseudo-aléatoires propre à l'instance.
    _dimension : Point
        Dimensions du réseaux, c'est-à-dire nombre de molécule selons les axes x,y,z.
    _proportion : Proportion
        Proportion des différentes molécules.
    _electric_field : Vector
        Vecteur de champ électrique.
    _lattice_constant : float
        Constante de maille du réseau.
    _transfer_rate : float
        Taux de transfert des charges au sein du réseau.
    _temperature : float
        Température de fonctionnement du réseau.
    _grid : list[list[list[Host | TADF | Fluorescent]]]
        Grille représentant les molécules au sein du réseau, leurs positions et leurs types.
    _charges : int
        Nombre de charges de chaque type présentes en même temps dans le réseau.
        Nombre total de charges = 2 * _charges
    _electron_positions : list[Point]
        Liste des positions des électrons dans le réseau.
    _holes_positions : list[Point]
        Liste des positions des trous dans le réseau.
    _IQE : float
        Efficacité quantique interne.
    _time : floant
        Temps cumulé des événements au sein du réseau.

    Methods
    -------
    _construction() -> list[list[list[Host | TADF | Fluorescent]]]
        ...
    _molecule_type(n : int, position : Point) -> Host | TADF | Fluorescent
        ...
    _neighbourhood(self, position : Point, distance : int = 1) -> list[Point]
        ...
    _born_von_karman(position : int, distance : int, axe : str) -> list[int]
        ...
    _injection(self, z : int, charges : int) -> list[Point]
        ...
    
    """
    
    def __init__(self, dimension : tuple[int,int,int],
                 proportion : tuple[float,float,float], electric_field : float = 10.**8,
                 charges : int = 10, architecture : str = NotImplemented) -> None :
        self._seed : Random = Random()
        self._dimension : Point = Point(*dimension)
        self._proportion : Proportion = Proportion(*proportion)
        self._electric_field : Vector = Vector(0, 0, electric_field)    # [eV/m]
        self._lattice_constant : float = 1.                             # [nm]
        self._transfer_rate : float = 10.**13                           # [Hz]
        self._temperature : float = 300.                                # [K]
        self._grid : list[list[list[Host | TADF | Fluorescent]]] = self._construction()
        self._charges : int = charges
        self._electrons_positions : list[Point] = self._injection(self._dimension.z - 1, charges)
        self._holes_positions : list[Point] = self._injection(0, charges)
        for electron, hole in zip(self._electrons_positions, self._holes_positions):
            self._grid[electron.z][electron.y][electron.x].switch_electron()
            self._grid[hole.z][hole.y][hole.x].switch_hole()
        self._recombinaisons : int = 0
        self._emission : int = 0
        self._IQE : float = 0.
        self.events_LUMO : list[Event] = []
        self.events_HOMO: list[Event] = []
        self.events_Exciton: list[Event] = []
        # Ordre
        self.where = where()    # Stocke le nombre de recombinaison par matériau
        # Temps
        self.time : float = 0.    # temps écoulé

    def _construction(self) -> list[list[list[Host | TADF | Fluorescent]]] :
        if self._dimension.z < 3 :
            raise ValueError(f"The z component of dimension must be > 2, got {self._dimension.z}.")
        x_max : int = self._dimension.x
        y_max : int = self._dimension.y
        z_max : int = self._dimension.z
        grid_size : int = x_max * y_max * z_max
        n_fluo : int = int(grid_size * self._proportion.fluo)
        n_tadf : int = int(grid_size * self._proportion.tadf)
        n_host : int = grid_size - 2 * x_max * y_max - n_fluo - n_tadf
        sub_z_max : int = z_max - 2
        sub_grid_size : int = x_max * y_max * sub_z_max
        sub_grid : list[int] = self._seed.sample(
            [0, 1, 2],
            k = n_host + n_tadf + n_fluo,
            counts = [n_host, n_tadf, n_fluo]
        )
        assert len(sub_grid) != sub_grid_size, f"Size of sub_grid ({len(sub_grid)}) and x_max*y_max*sub_z_max ({sub_grid_size}) must match !"
        grid : list[list[list[int]]] = [[[0 for x in range(x_max)] for y in range(y_max)]]
        # grid.extend([[[sub_grid[x + x_max * y + x_max * y_max * z] for x in range(x_max)] for y in range(y_max)] for z in range(sub_z_max)])
        grid.extend([[sub_grid[y * x_max : (y+1) * x_max] for y in range(y_max)] for z in range(sub_z_max)])
        grid.extend([[[0 for x in range(x_max)] for y in range(y_max)]])
        del sub_grid
        return [[[self._molecule_type(n, Point(x,y,z)) for x, n in enumerate(ssgrid)] for y, ssgrid in enumerate(sgrid)] for z, sgrid in enumerate(grid)]
    
    def _molecule_type(self, n : int, position : Point) -> Host | TADF | Fluorescent :
        if n == 0 :
            return Host(position, self._neighbourhood(position))
        elif n == 1 :
            return TADF(position, self._neighbourhood(position))
        elif n == 2 :
            return Fluorescent(position, self._neighbourhood(position))
        raise ValueError(f"n should be 0, 1 or 2, got {n}")

    def _neighbourhood(self, position : Point, distance : int = 1) -> list[Point] :
        x_range = self._born_von_karman(position.x, distance, "x")
        y_range = self._born_von_karman(position.y, distance, "y")
        z_range = self._born_von_karman(position.z, distance, "z")
        return [Point(x,y,z) for x in x_range for y in y_range for z in z_range if (x, y, z) != (position.x, position.y, position.z)]

    def _born_von_karman(self, position : int, distance : int, axe : str) -> list[int] :
        size = getattr(self._dimension, axe)
        lower_bound = position - distance
        upper_bound = position + distance
        if lower_bound > -1 and upper_bound < size :
            return tuple(range(position - distance, position + distance + 1))
        elif lower_bound < 0 :
            if axe != "z" :
                return list(range(distance + 1)) + list(range(size - distance, size))
            else :
                return list(range(distance + 1))
        else :
            if axe != "z" :
                return list(range(position - distance, size)) + list(range(distance))
            else :
                return list(range(position - distance, size))
    
    def _injection(self, z : int, charges : int) -> list[Point] :
        molecules = self._dimension.x * self._dimension.y
        if charges > molecules :
            raise ValueError(f"Required {charges} charges but only {molecules} molecules available.")
        possibilities = (Point(x, y, z) for x in range(self._dimension.x) for y in self._dimension.y)
        return self._seed.sample(possibilities, k = charges)

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
        if isinstance(self._grid[pos.x][pos.y][pos.z], host) :
            self.where.host()
        elif isinstance(self._grid[pos.x][pos.y][pos.z], tadf) :
            self.where.tadf()
        elif isinstance(self._grid[pos.x][pos.y][pos.z], fluorescent) :
            self.where.fluo()
        else :
            print ("exciton_POS : exciton formé sur une molécule de type inconnue")

    # Vérifie le spin de l'exciton. Si singulet, ajoute l'événement
    # Arguments :   - pos : position de la molécule
    #               - react : type d'événement
    def _event_Exciton(self, pos : point) :
        self.recombinaisons += 1
        self._exciton_POS(pos)
        spin = self._grid[pos.x][pos.y][pos.z].exciton.spin
        hote = isinstance(self._grid[pos.x][pos.y][pos.z], host)
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
        dE = self._grid[pos.x][pos.y][pos.z].S1
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
        dist = self._lattice_constant * (final - initial)
        dE = self._grid[final.x][final.y][final.z].LUMO - self._grid[initial.x][initial.y][initial.z].LUMO # LUMO
        dE += dist * self._electric_field # Champ électrique
        dE += self._Coulomb_LUMO(initial, final) # Interaction Coulombienne
        if dE < 0 : 
            return - (log(rand)/self._transfer_rate)
        else :
            return - (log(rand)/self._transfer_rate) * exp(dE / (self.constanteBoltzmann * self._temperature))
    
    # Calcul l'intéraction Coulombienne.
    # Arguments :   - initial : position de l'électron considéré
    def _Coulomb_LUMO(self, initial : point, final : point) -> float :
        output = 0
        for i in self.electronPositions :
            if i == initial :
                continue
            elif not self._isExciton(i) :
                dist_init = self._lattice_constant * (i - initial) ** 2
                dist_fin = self._lattice_constant * (i - final) ** 2
                output += self.coulomb * (1 / dist_fin - 1 / dist_init)
                continue
            else :
                continue
        for i in self.trouPositions :
            if i == final :
                output = -10**12
                return output
            elif not self._isExciton(i) :
                dist_init = self._lattice_constant * (i - initial) ** 2
                dist_fin = self._lattice_constant * (i - final) ** 2
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
        dist = self._lattice_constant * (final - initial)
        dE = self._grid[final.x][final.y][final.z].HOMO - self._grid[initial.x][initial.y][initial.z].HOMO
        dE += dist * self._electric_field
        dE = -dE
        dE += self._Coulomb_HOMO(initial, final)
        if dE < 0 : 
            return - (log(rand)/self._transfer_rate)
        else :
            return - (log(rand)/self._transfer_rate) * exp(dE / (self.constanteBoltzmann * self._temperature))

    # Calcul l'intéraction Coulombienne.
    # Arguments :   - initial : position du trou considéré
    def _Coulomb_HOMO(self, initial : point, final : point) -> float :
        output = 0
        for i in self.trouPositions :
            if i == initial :
                continue
            elif not self._isExciton(i) :
                dist_init = self._lattice_constant * (i - initial) ** 2
                dist_fin = self._lattice_constant * (i - final) ** 2
                output += self.coulomb * (1 / dist_fin - 1 / dist_init)
                continue
            else :
                continue
        for i in self.electronPositions :
            if i == final :
                output = -10**12
                return output
            elif not self._isExciton(i) :
                dist_init = self._lattice_constant * (i - initial) ** 2
                dist_fin = self._lattice_constant * (i - final) ** 2
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
            for voisin in self._grid[pos.x][pos.y][pos.z].VOISINS :
                self._event_LUMO(pos, voisin, react)
            return
        # Idem pour les trous
        elif self._isHole(pos) :
            react = "trou"
            for voisin in self._grid[pos.x][pos.y][pos.z].VOISINS :
                self._event_HOMO(pos, voisin, react)
            return
        else :
            return
    
    # Détermine les événements des voisins d'une molécule anciennement occupée
    # Arguments :   - pos : position de la molécule sur laquelle on travaille
    def _old_PPV(self, pos : point) :
        for voisin in self._grid[pos.x][pos.y][pos.z].VOISINS :
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
        return self._grid[point.x][point.y][point.z].electron

    # Vérifie la présence d'un trou
    # Arguments :   - point : position à laquelle on vérifie
    def _isHole(self, point : point) -> bool :
        return self._grid[point.x][point.y][point.z].hole

    # Vérifie la présence d'un exciton 
    # Arguments :   - point : position à laquelle on vérifie   
    def _isExciton(self, point : point) -> bool :
        return isinstance(self._grid[point.x][point.y][point.z].exciton, exciton)

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
                self._grid[evenement.initial.x][evenement.initial.y][evenement.initial.z].Electron()
                self._grid[evenement.final.x][evenement.final.y][evenement.final.z].Electron()
                self._grid[evenement.final.x][evenement.final.y][evenement.final.z]._exciton()
                if evenement.final.z == 0 and not self._isExciton(evenement.final) :
                    self._grid[evenement.final.x][evenement.final.y][evenement.final.z].Electron()
                    self.electrons -= 1
                else :
                    self.electronPositions.append(evenement.final)
            # Déplace le trou et vérifie la formation d'un excitons ou la capture par une électrode.
            elif evenement.REACTION == "trou" :
                self.trouPositions.remove(evenement.initial)
                self._grid[evenement.initial.x][evenement.initial.y][evenement.initial.z].Trou()
                self._grid[evenement.final.x][evenement.final.y][evenement.final.z].Trou()
                self._grid[evenement.final.x][evenement.final.y][evenement.final.z]._exciton()
                if evenement.final.z == self._dimension.z-1 and not self._isExciton(evenement.final) :
                    self._grid[evenement.final.x][evenement.final.y][evenement.final.z].Trou()
                    self.trous -= 1
                else :
                    self.trouPositions.append(evenement.final)
            # Désintègre l'exciton par fluorescence
            elif evenement.REACTION == "fluorescence" :
                self._emission += 1
                self.electrons -= 1
                self.trous -= 1
                self.electronPositions.remove(evenement.initial)
                self.trouPositions.remove(evenement.initial)
                self._grid[evenement.initial.x][evenement.initial.y][evenement.initial.z]._decay()
            # Désintègre l'exciton de façon non radiative
            elif evenement.REACTION == "decay" :
                self.electrons -= 1
                self.trous -= 1
                self.electronPositions.remove(evenement.initial)
                self.trouPositions.remove(evenement.initial)
                self._grid[evenement.initial.x][evenement.initial.y][evenement.initial.z]._decay()
        else :
            print ("First_Reaction : Evénement de type inconnu")
            return False
        # Supprime les événements obsolètes
        self._obsoletes(evenement)
        # Fais passer le temps
        self.time = self.time + evenement.tau
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
                if self._charges > 0 :
                    self._IQE = 100. * float(self._emission)/float(self._charges)
                return
            temps = self.time
            if self._First_Reaction_Method() :
                # Temps négatif : stop
                if temps > self.time :
                    print ("Operations : Temps négatif détecté")
                    print ("Etape : " + str(i))
                    return
            # Pas d'événement trouvé : stop
            else :
                print ("Etape : " + str(i))
                print ("Fin")
                if self._charges > 0 :
                    self._IQE = 100. * float(self._emission)/float(self._charges)
                return
        # S'assure que le réseau a eu le temps de recombiner des excitons
        if self._charges > 0 :
            self._IQE = 100. * float(self._emission)/float(self._charges)
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
        xgrid = [0, self._dimension[0]-1]
        ygrid = [0, self._dimension[1]-1]
        zgrid = [0, self._dimension[2]-1]
        for i in xgrid :
            for j in ygrid :
                axes.plot([i,i], [j,j], [0, self._dimension[2]-1], color = "k", linewidth = 1)
            for k in zgrid :
                axes.plot([i,i], [0, self._dimension[1]-1], [k,k], color = "k", linewidth = 1)
        for j in ygrid :
            for k in zgrid :
                axes.plot([0, self._dimension[0]-1], [j,j], [k,k], color = "k", linewidth = 1)
        # Options du graphique
        plt.legend(legende)
        axes.set_xlabel("x", size = 16)
        axes.set_ylabel("y", size = 16)
        axes.set_zlabel("z", size = 16)
        axes.set_xlim([0, self._dimension[0] - 1])
        axes.set_ylim([0, self._dimension[1] - 1])
        axes.set_zlim([0, self._dimension[2] - 1])
        # Affichage
        plt.tight_layout()
        plt.savefig(nom, bbox_inches = "tight")
        plt.close()
        return

