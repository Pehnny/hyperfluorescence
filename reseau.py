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
    
    ###############################################
    ####____Méthodes de démarrage du réseau____####
    ###############################################
    def __init__(self, dimension : tuple[int,int,int],
                 proportions : tuple[float,float,float], electric_field : float = 10.**8,
                 charges : int = 10, architecture : str = NotImplemented) -> None :
        self._seed : Random = Random()
        self._lattice_parameters_creation(dimension, proportions, electric_field, charges)
        self._grid : list[list[list[Host | TADF | Fluorescent]]] = self._construction()
        self._charges_injection()
        self._events_creation()
        self._recombinaisons : int = 0
        self._emission : int = 0
        self._IQE : float = 0.
        self.time : float = 0.

    def _lattice_parameters_creation(self, dimension : tuple[int,int,int],
                                     proportions : tuple[float,float,float], electric_field : float,
                                     charges : int) -> None :
        self._dimension : Point = Point(*dimension)
        self._proportions : Proportion = Proportion(*proportions)
        self._electric_field : Vector = Vector(0, 0, electric_field)    # [eV/m]
        self._lattice_constant : float = 1.                             # [nm]
        self._transfer_rate : float = 10.**13                           # [Hz]
        self._temperature : float = 300.                                # [K]
        self._charges : int = charges

    def _construction(self) -> list[list[list[Host | TADF | Fluorescent]]] :
        """Méthode construisant le réseau. Les couches en z = 0 et z = z_max sont entièrement composées
        de molécules hôte. Les autres couches sont composées aléatoirement selon les proportions de chaque
        type de molécule.
        """
        if self._dimension.z < 3 :
            raise ValueError(f"The z component of dimension must be > 2, got {self._dimension.z}.")
        x_max : int = self._dimension.x
        y_max : int = self._dimension.y
        z_max : int = self._dimension.z
        grid_size : int = x_max * y_max * z_max
        n_fluo : int = int(grid_size * self._proportions.fluo)
        n_tadf : int = int(grid_size * self._proportions.tadf)
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
    
    def _molecule_type(self, n : int,
                       position : Point) -> Host | TADF | Fluorescent :
        """Méthode convertissant les entiers en type de molécule.
        """
        if n == 0 :
            return Host(position, self._neighbourhood(position))
        elif n == 1 :
            return TADF(position, self._neighbourhood(position))
        elif n == 2 :
            return Fluorescent(position, self._neighbourhood(position))
        raise ValueError(f"n should be 0, 1 or 2, got {n}")

    def _neighbourhood(self, position : Point,
                       distance : int = 1) -> list[Point] :
        """Méthode déterminant les voisins d'une molécule.
        """
        x_range = self._born_von_karman(position.x, distance, "x")
        y_range = self._born_von_karman(position.y, distance, "y")
        z_range = self._born_von_karman(position.z, distance, "z")
        return [Point(x,y,z) for x in x_range for y in y_range for z in z_range if (x, y, z) != (position.x, position.y, position.z)]

    def _born_von_karman(self, position : int,
                         distance : int, axe : str) -> list[int] :
        """Méthode déterminant quels voisins sont valides selon le principe de motif périodique.
        L'axe z est le seul qui ne reboucle pas.
        """
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
    
    def _charges_injection(self) -> None :
        molecules = self._dimension.x * self._dimension.y
        if self._charges > molecules :
            raise ValueError(f"Required {self._charges} charges but only {molecules} molecules available.")
        positions = (
                Point(x, y, 0)
                for x in range(self._dimension.x)
                for y in range(self._dimension.y)
            )
        self._electrons_positions : list[Point] = self._seed.sample(positions, k = self._charges)
        positions = (
                Point(x, y, self._dimension.z - 1)
                for x in range(self._dimension.x)
                for y in range(self._dimension.y)
            )
        self._holes_positions : list[Point] = self._seed.sample(positions, k = self._charges)
        self._exciton_positions : list[Point] = []
        for electron, hole in zip(self._electrons_positions, self._holes_positions) :
            self._grid[electron.z][electron.y][electron.x].switch_electron()
            self._grid[hole.z][hole.y][hole.x].switch_hole()
    
    def _events_creation(self) -> None :
        self._move_electron_events : list[Event] = self._init_move_electron()
        self._move_hole_events : list[Event] = self._init_move_hole()
        self._move_exciton_events : list[Event] = []
        self._decay_events : list[Event] = []
        self._isc_events : list[Event] = []
        self._capture_events : list[Event] = []
        self._unbound_events : list[Event] = []

    def _init_move_electron(self) -> list[Event] :
        output : list[Event] = []
        for position in self._electrons_positions :
            molecule = self.get_molecule(position)
            output.extend((
                Event(position, neighbour, self._LUMO(position, neighbour), EVENTS["move"]["electron"]) 
                for neighbour in molecule.neighbors
                if not self.get_molecule(neighbour).electron    
            ))
        return output
    
    def _init_move_hole(self) -> list[Event] :
        output : list[Event] = []
        for position in self._holes_positions :
            molecule = self.get_molecule(position)
            output.extend((
                Event(position, neighbour, self._LUMO(position, neighbour), EVENTS["move"]["electron"]) 
                for neighbour in molecule.neighbors
                if not self.get_molecule(neighbour).hole  
            ))
        return output
    

    #########################################################
    ####____Méthodes de tri des événements inutulisés____####
    #########################################################
    def _remove_unused_event(self, event : Event) -> None :
        if is_move_event(event.kind) :
            if is_electron(event.particule) :
                self._move_electron_events.remove(event)
            elif is_hole(event.particule) :
                self._move_hole_events.remove(event)
            elif is_exciton(event.particule) :
                self._move_exciton_events.remove(event)
            raise ValueError(f"event.particule is {event.particule}, expected one of {PARTICULES}.")
        elif is_decay_event(event.kind) :
            self._decay_events.remove(event)
        elif is_ISC_event(event.kind) :
            ...
        elif is_unbound_event(event.kind) :
            ...
        elif is_capture_event(event.kind) :
            self._capture_events.remove(event)
        raise TypeError(f"event must be of type Event, got {type(event)}.")

    ####____Méthode d'extraction des données____####
    # def _exciton_POS(self, pos : point) :
    #     if isinstance(self._grid[pos.x][pos.y][pos.z], host) :
    #         self.where.host()
    #     elif isinstance(self._grid[pos.x][pos.y][pos.z], tadf) :
    #         self.where.tadf()
    #     elif isinstance(self._grid[pos.x][pos.y][pos.z], fluorescent) :
    #         self.where.fluo()
    #     else :
    #         print ("exciton_POS : exciton formé sur une molécule de type inconnue")

    ####################################################
    ####____Méthodes qui génèrent les événements____####
    ####################################################
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
    
    def get_molecule(self, position : Point) -> Host | TADF | Fluorescent :
        return self._grid[position.z][position.y][position.x]

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
    
    ##################################################################
    ####____Méthode d'affichage du réseau dans son état actuel____####
    ##################################################################
    def get_molecule_type(self, position : Point) -> type :
        """Méthode retournant le type de la molécule à la position donnée.
        """
        return type(self._grid[position.z][position.y][position.x])


    def Plot(self, nom : str = "Inconnu") :
        use("Agg")
        plt.figure(dpi=50)
        axes = plt.axes(projection = "3d")

        color = "b"
        electron_host = [
            position
            for position in self._electrons_positions
            if self.get_molecule_type(position) is Host
        ]
        if len(electron_host) > 0 :
            marker_style = "o"
            x = (position.x for position in electron_host)
            y = (position.y for position in electron_host)
            z = (position.z for position in electron_host)
            axes.scatter(x, y, z, c = color, marker = marker_style)
        electron_tadf = [
            position
            for position in self._electrons_positions
            if self.get_molecule_type(position) is TADF
        ]
        if len(electron_tadf) > 0 :
            marker_style = "s"
            x = (position.x for position in electron_tadf)
            y = (position.y for position in electron_tadf)
            z = (position.z for position in electron_tadf)
            axes.scatter(x, y, z, c = color, marker = marker_style)
        electron_fluorescent = [
            position
            for position in self._electrons_positions
            if self.get_molecule_type(position) is Fluorescent
        ]
        if len(electron_fluorescent) > 0 :
            marker_style = "^"
            x = (position.x for position in electron_fluorescent)
            y = (position.y for position in electron_fluorescent)
            z = (position.z for position in electron_fluorescent)
            axes.scatter(x, y, z, c = color, marker = marker_style)

        color = "r"
        hole_host = [
            position
            for position in self._holes_positions
            if self.get_molecule_type(position) is Host
        ]
        if len(hole_host) > 0 :
            marker_style = "o"
            x = (position.x for position in hole_host)
            y = (position.y for position in hole_host)
            z = (position.z for position in hole_host)
            axes.scatter(x, y, z, c = color, marker = marker_style)
        hole_tadf = [
            position
            for position in self._holes_positions
            if self.get_molecule_type(position) is TADF
        ]
        if len(hole_tadf) > 0 :
            marker_style = "s"
            x = (position.x for position in hole_host)
            y = (position.y for position in hole_host)
            z = (position.z for position in hole_host)
            axes.scatter(x, y, z, c = color, marker = marker_style)
        hole_fluorescent = [
            position
            for position in self._holes_positions
            if self.get_molecule_type(position) is Fluorescent
        ]
        if len(hole_fluorescent) > 0 :
            marker_style = "^"
            x = (position.x for position in hole_host)
            y = (position.y for position in hole_host)
            z = (position.z for position in hole_host)
            axes.scatter(x, y, z, c = color, marker = marker_style)

        x_grid = [0, self._dimension[0]-1]
        y_grid = [0, self._dimension[1]-1]
        z_grid = [0, self._dimension[2]-1]
        for x in x_grid :
            for y in y_grid :
                axes.plot([x,x], [y,y], [0, self._dimension[2]-1], linestyle = "dashed", color = "k")
            for z in z_grid :
                axes.plot([x,x], [0, self._dimension[1]-1], [z,z], linestyle = "solid", color = "k")
        for y in y_grid :
            for z in z_grid :
                axes.plot([0, self._dimension[0]-1], [y,y], [z,z], linestyle = "solid", color = "k")

        axes.set_xlabel("x", size = 16)
        axes.set_ylabel("y", size = 16)
        axes.set_zlabel("z", size = 16)
        axes.set_xlim([0, self._dimension.x - 1])
        axes.set_ylim([0, self._dimension.y - 1])
        axes.set_zlim([0, self._dimension.z - 1])
        plt.tight_layout()
        plt.savefig(nom, bbox_inches = "tight")
        plt.close()
        return

