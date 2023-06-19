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

from math import exp, log, prod
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


class Lattice :
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
    _holes_locations : list[Point]
        Liste des positions des trous dans le réseau.
    _IQE : float
        Efficacité quantique interne.
    _time : floant
        Temps cumulé des événements au sein du réseau.

    Methods
    -------
    _lattice_creation() -> list[list[list[Host | TADF | Fluorescent]]]
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
    def __init__(self, dimension : tuple[int,int,int], proportions : tuple[float,float,float],
                 electric_field : float = 10.**(-1), charges : int = 10, architecture : str = NotImplemented) -> None :
        self._seed : Random = Random()
        self._lattice_parameters_creation(dimension, proportions, electric_field, charges)
        self._grid : list[list[list[Host | TADF | Fluorescent]]] = self._lattice_creation()
        self._charges_injection()
        self._events_creation()
        self._recombinations : int = 0
        self._emission : int = 0
        self._IQE : float = 0.
        self._step : int = 0
        self.time : float = 0.

    def _lattice_parameters_creation(self, dimension : tuple[int,int,int], proportions : tuple[float,float,float],
                                     electric_field : float, charges : int) -> None :
        if len(dimension) != 3 :
            raise IndexError(f"dimensions length must be 3, got {len(dimension)}")
        if dimension[-1] < 3 :
            raise ValueError(f"The last component of dimension must be > 2, got {dimension[-1]}.")
        if len(proportions) != 3 :
            raise IndexError(f"proportions length must be 3, got {len(dimension)}")
        if sum(proportions) != 1. :
            norm = sum(dimension)
            proportions = (dimension[0] / norm, dimension[1] / norm, dimension[2] / norm)
        if  min(proportions) * prod(dimension) < 1 :
            minimum = 1 / min(proportions)
            raise ValueError(f"Not enough molecules for current proportions. Expected at least {minimum}, got {prod(dimension)}")
        self._dimension : Point = Point(*dimension)
        self._proportions : Proportion = Proportion(*proportions)
        self._electric_field : Vector = Vector(0, 0, electric_field)    # [eV/nm]
        self._lattice_constant : float = 1.                             # [nm]
        self._transfer_rate : float = 10.**13                           # [Hz]
        self._temperature : float = 300.                                # [K]
        self._charges : int = charges

    def _lattice_creation(self) -> list[list[list[Host | TADF | Fluorescent]]] :
        x_max : int = self._dimension.x
        y_max : int = self._dimension.y
        z_max : int = self._dimension.z
        grid_size : int = x_max * y_max * z_max
        n_fluo : int = max(int(grid_size * self._proportions.fluo), 1)
        n_tadf : int = max(int(grid_size * self._proportions.tadf), 1)
        n_host : int = max(int(grid_size - 2 * x_max * y_max - n_fluo - n_tadf), 1)
        sub_z_max : int = z_max - 2
        sub_grid_size : int = x_max * y_max * sub_z_max
        sub_grid : list[int] = self._seed.sample(
            [0, 1, 2],
            k = n_host + n_tadf + n_fluo,
            counts = [n_host, n_tadf, n_fluo]
        )
        assert len(sub_grid) == sub_grid_size, f"Size of sub_grid ({len(sub_grid)}) and x_max*y_max*sub_z_max ({sub_grid_size}) must match !"
        grid : list[list[list[int]]] = [[[0 for x in range(x_max)] for y in range(y_max)]]
        grid.extend([[sub_grid[y * x_max : (y+1) * x_max] for y in range(y_max)] for z in range(sub_z_max)])
        grid.extend([[[0 for x in range(x_max)] for y in range(y_max)]])
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
        size : int = getattr(self._dimension, axe)
        lower_bound : int = position - distance
        upper_bound : int = position + distance
        if lower_bound > -1 and upper_bound < size :
            return list(range(position - distance, position + distance + 1))
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
        self._electrons_locations : list[Point] = []
        self._holes_locations : list[Point] = []
        self._excitons_locations : list[Point] = []
        molecules : int = self._dimension.x * self._dimension.y
        if self._charges > molecules :
            raise ValueError(f"Required {self._charges} charges but only {molecules} molecules available.")
        positions = [
                Point(x, y, self._dimension.z - 1)
                for x in range(self._dimension.x)
                for y in range(self._dimension.y)
            ]
        self._electrons_locations.extend(self._seed.sample(positions, k = self._charges))
        positions = [
                Point(x, y, 0)
                for x in range(self._dimension.x)
                for y in range(self._dimension.y)
            ]
        self._holes_locations.extend(self._seed.sample(positions, k = self._charges))
        for electron, hole in zip(self._electrons_locations, self._holes_locations) :
            self._grid[electron.z][electron.y][electron.x].switch_electron()
            self._grid[hole.z][hole.y][hole.x].switch_hole()
        
    
    def _events_creation(self) -> None :
        self._move_electron_events : list[Event] = self._init_move_electron_events()
        self._move_hole_events : list[Event] = self._init_move_hole_events()
        self._move_exciton_events : list[Event] = []
        self._decay_events : list[Event] = []
        self._isc_events : list[Event] = []
        self._binding_events : list[Event] = []
        self._capture_events : list[Event] = []

    def _init_move_electron_events(self) -> list[Event] :
        molecules = (
            self._get_molecule(position)
            for position in self._electrons_locations
        )
        events = (
            Event(molecule.position, neighbour, self._time_move_electron(molecule.position, neighbour), EVENTS["move"], PARTICULES["electron"])
            for molecule in molecules
            for neighbour in molecule.neighbourhood
            if neighbour not in self._electrons_locations
        )
        return list(events)

    def _init_move_hole_events(self) -> list[Event] :
        molecules = (
            self._get_molecule(position)
            for position in self._holes_locations
        )
        events = (
            Event(molecule.position, neighbour, self._time_move_hole(molecule.position, neighbour), EVENTS["move"], PARTICULES["hole"])
            for molecule in molecules
            for neighbour in molecule.neighbourhood
            if neighbour not in self._holes_locations
        )
        return list(events)



    ##################################################################
    ####____Méthodes de calcul des durées de chaque événement_____####
    ##################################################################
    def _time_move_electron(self, initial : Point, final : Point) -> float :
        rng = 1. - self._seed.random()
        movement : Vector = (final - initial) * self._lattice_constant
        delta_energy : float = self._lumo_energy(initial, final)
        delta_energy += self._electric_field * movement
        delta_energy += self._electron_electrostatic_energy(initial, final)
        if delta_energy < 0 : 
            return - log(rng) / self._transfer_rate
        else :
            return - (log(rng) / self._transfer_rate) * exp(delta_energy / (cst.BOLTZMANN * self._temperature))
        
    def _lumo_energy(self, initial : Point, final : Point) -> float :
        return self._get_molecule(final).lumo_energy - self._get_molecule(initial).lumo_energy
    
    def _electron_electrostatic_energy(self, initial : Point, final : Point) -> float :
        output : float = 0.
        for location in self._electrons_locations :
            if location == initial : continue
            old_r : Vector = (location - initial) * self._lattice_constant
            new_r : Vector = (location - final) * self._lattice_constant
            delta_r : float = 1. / new_r.norm() - 1. / old_r.norm()
            output += cst.ELECTROSTATIC * delta_r
        for location in self._holes_locations :
            if location == final : continue
            old_r = (location - initial) * self._lattice_constant
            new_r = (location - final) * self._lattice_constant
            delta_r = 1. / new_r.norm() - 1. / old_r.norm()
            output -= cst.ELECTROSTATIC * delta_r
        return output

    def _time_move_hole(self, initial : Point, final : Point) -> float :
        rng = 1. - self._seed.random()
        movement : Vector = (final - initial) * self._lattice_constant
        delta_energy = self._homo_energy(initial, final)
        delta_energy -= self._electric_field * movement
        delta_energy += self._hole_electrostatic_energy(initial, final)
        if delta_energy < 0 : 
            return - log(rng) / self._transfer_rate
        else :
            return - (log(rng) / self._transfer_rate) * exp(delta_energy / (cst.BOLTZMANN * self._temperature))
        
    def _homo_energy(self, initial : Point, final : Point) -> float :
        return self._get_molecule(initial).homo_energy - self._get_molecule(final).homo_energy

    # Calcul l'intéraction Coulombienne.
    # Arguments :   - initial : position du trou considéré
    def _hole_electrostatic_energy(self, initial : Point, final : Point) -> float :
        output : float = 0.
        for location in self._holes_locations :
            if location == initial : continue
            old_r : Vector = (location - initial) * self._lattice_constant
            new_r : Vector = (location - final) * self._lattice_constant
            delta_r : float = 1. / new_r.norm() - 1. / old_r.norm()
            output += cst.ELECTROSTATIC * delta_r
        for location in self._electrons_locations :
            if location == final : continue
            old_r = (location - initial) * self._lattice_constant
            new_r = (location - final) * self._lattice_constant
            delta_r = 1. / new_r.norm() - 1. / old_r.norm()
            output -= cst.ELECTROSTATIC * delta_r
        return output

    

    ################################################################################
    ####____Méthodes qui suppriment les événements qui ne sont plus utilisés____####
    ################################################################################
    def _remove_move_electron_events(self, event : Event) -> None :
        while True :
            try :
                self._move_electron_events.remove(event)
            except ValueError :
                break
        
    def _remove_move_hole_events(self, event : Event) -> None :
        while True :
            try :
                self._move_hole_events.remove(event)
            except ValueError :
                break

    def _remove_bound_event(self, event : Event) -> None :
        self._binding_events.remove(event)

    def _remove_decay_event(self, event : Event) -> None :
        self._decay_events.remove(event)

    def _remove_capture_event(self, event : Event) -> None :
        self._capture_events.remove(event)



    ############################################################################
    ####____Méthodes qui génèrent les nouveaux événements à chaque étape____####
    ############################################################################
    def _new_move_electron_events(self, position : Point) -> None :
        neighbourhood = self._get_molecule(position).neighbourhood
        events = (
            Event(position, neighbour, self._time_move_electron(position, neighbour), EVENTS["move"], PARTICULES["electron"])
            for neighbour in neighbourhood
            if not self._get_molecule(neighbour).electron
        )
        self._move_electron_events.extend(events)

    def _new_move_hole_events(self, position : Point) -> None :
        neighbourhood = self._get_molecule(position).neighbourhood
        events = (
            Event(position, neighbour, self._time_move_hole(position, neighbour), EVENTS["move"], PARTICULES["hole"])
            for neighbour in neighbourhood
            if not self._get_molecule(neighbour).hole
        )
        self._move_hole_events.extend(events)

    def _new_bound_event(self, position : Point) -> None :
        event = Event(position, position, 0., EVENTS["bound"], PARTICULES["exciton"])
        self._binding_events.append(event)

    def _new_decay_event(self, position : Point) -> None :
        event = Event(position, position, 0., EVENTS["decay"], PARTICULES["exciton"])
        self._decay_events.append(event)

    def _new_capture_electron_event(self, position : Point) -> None :
        event = Event(position, position, 0., EVENTS["capture"], PARTICULES["electron"])
        self._capture_events.append(event)

    def _new_capture_hole_event(self, position : Point) -> None :
        event = Event(position, position, 0., EVENTS["capture"], PARTICULES["hole"])
        self._capture_events.append(event)



    ####################################################
    ####____Méthodes de transformation du réseau____####
    ####################################################
    def _move_electron(self, initial : Point, final : Point) -> None :
        self._grid[initial.z][initial.y][initial.x].switch_electron()
        self._grid[final.z][final.y][final.x].switch_electron()
        self._electrons_locations.remove(initial)
        self._electrons_locations.append(final)

    def _move_hole(self, initial : Point, final : Point) -> None :
        self._grid[initial.z][initial.y][initial.x].switch_hole()
        self._grid[final.z][final.y][final.x].switch_hole()
        self._holes_locations.remove(initial)
        self._holes_locations.append(final)
    
    def _form_exciton(self, position : Point) -> None :
        self._grid[position.z][position.y][position.x].generate_exciton()
        self._electrons_locations.remove(position)
        self._holes_locations.remove(position)
        self._excitons_locations.append(position)

    def _capture_electron(self, position : Point) -> None :
        self._grid[position.z][position.y][position.x].switch_electron()
        self._electrons_locations.remove(position)

    def _capture_hole(self, position : Point) -> None :
        self._grid[position.z][position.y][position.x].switch_hole()
        self._holes_locations.remove(position)

    def _electron_reinjection(self) -> None :
        positions = [
                Point(x, y, self._dimension.z - 1)
                for x in range(self._dimension.x)
                for y in range(self._dimension.y)
                if Point(x, y, self._dimension.z - 1) not in self._electrons_locations
            ]
        self._electrons_locations.append(self._seed.choice(positions))
        position = self._electrons_locations[-1]
        self._grid[position.z][position.y][position.x].switch_electron()

    def _hole_reinjection(self) -> None :
        positions = [
                Point(x, y, 0)
                for x in range(self._dimension.x)
                for y in range(self._dimension.y)
                if Point(x, y, 0) not in self._holes_locations
            ]
        self._holes_locations.append(self._seed.choice(positions))
        position = self._holes_locations[-1]
        self._grid[position.z][position.y][position.x].switch_hole()

    def _decay(self, position : Point) -> None :
        photon = self._grid[position.z][position.y][position.x].exciton_decay()
        self._excitons_locations.remove(position)
        self._recombinations += 1
        if photon : self._emission += 1



    ##################################################
    ####____Algorithme "First Reaction Method"____####
    ##################################################
    def _first_reaction_method(self) -> bool :
        #   Vérifie si il reste un événement et récupère le plus rapide
        events : list[Event] = self._get_all_events()
        events.sort(reverse = True)
        try : event = events.pop()
        except IndexError : return False
        #   Traite les événements de type "move"
        if event.kind == EVENTS["move"] :
            #   Traite le cas d'un électron
            if event.particule == PARTICULES["electron"] :
                self._remove_move_electron_events(event)
                self._move_electron(event.initial, event.final)
                molecule = self._get_molecule(event.final)
                if not molecule.hole and event.final.z != 0 :
                    self._update_events(event.tau)
                    self._new_move_electron_events(event.final)
                elif molecule.hole :
                    event = Event(event.final, event.final, 0., EVENTS["move"], PARTICULES["hole"])
                    self._remove_move_hole_events(event)
                    self._update_events(event.tau)
                    self._new_bound_event(event.final)
                elif event.final.z == 0 :
                    self._update_events(event.tau)
                    self._new_capture_electron_event(event.final)
            #   Traite le cas d'un trou
            elif event.particule == PARTICULES["hole"] :
                self._remove_move_hole_events(event)
                self._move_hole(event.initial, event.final)
                molecule = self._get_molecule(event.final)
                if not molecule.electron and event.final.z != (self._dimension.z - 1) :
                    self._update_events(event.tau)
                    self._new_move_hole_events(event.final)
                elif molecule.electron :
                    event = Event(event.final, event.final, 0., EVENTS["move"], PARTICULES["electron"])
                    self._remove_move_electron_events(event)
                    self._update_events(event.tau)
                    self._new_bound_event(event.final)
                elif event.final.z == (self._dimension.z - 1) :
                    self._update_events(event.tau)
                    self._new_capture_hole_event(event.final)
            #   Traite le cas d'un exciton (non implémenté)
            elif event.particule == PARTICULES["exciton"] :
                ...
        #   Traite les événements de type formation d'exciton
        elif event.kind == EVENTS["bound"] :
            self._remove_bound_event(event)
            self._form_exciton(event.final)
            self._new_decay_event(event.final)
            ... # move ou unbound si Host ; ISC ou Forster si TADF ; Decay si Fluorescent
        #   Traite les événements de type conversion intersystème
        elif event.kind == EVENTS["ISC"] :
            ...
        #   Traite les événements de type transfert d'énergie de Forster
        elif event.kind == EVENTS["Forster"] :
            ...
        #   Traite les événements de type recombinaison
        elif event.kind == EVENTS["decay"] :
            self._remove_decay_event(event)
            self._decay(event.initial)
            self._electron_reinjection()
            self._new_move_electron_events(self._electrons_locations[-1])
            self._hole_reinjection()
            self._new_move_hole_events(self._holes_locations[-1])
        #   Traite les événements de type séparation d'exciton :
        elif event.kind == EVENTS["unbound"] :
            ...
        elif event.kind == EVENTS["capture"] :
            if event.particule == PARTICULES["electron"] :
                self._remove_capture_event(event)
                self._capture_electron(event.final)
                self._electron_reinjection()
                self._new_move_electron_events(self._electrons_locations[-1])
            elif event.particule == PARTICULES["hole"] :
                self._remove_capture_event(event)
                self._capture_hole(event.final)
                self._hole_reinjection()
                self._new_move_hole_events(self._holes_locations[-1])
        return True
    
    def _update_events(self, time : float) -> None :
        for event in self._move_electron_events :
            event.tau -= time
        for event in self._move_hole_events :
            event.tau -= time
        for event in self._move_exciton_events :
            event.tau -= time
        for event in self._binding_events :
            event.tau -= time
        for event in self._decay_events :
            event.tau -= time
        for event in self._isc_events :
            event.tau -= time
        for event in self._capture_events :
            event.tau -= time

    def operations(self, stop : int, start : int = 0) -> None :
        for i in range(start, stop) :
            #   Affiche l'étape régulièrement
            if not self._step%1000 :
                print (f"Step {self._step}")
            self._step += 1
            time = self.time
            #   Exécute l'évenement suivant et s'assure que le temps n'a pas diminué.
            if self._first_reaction_method() :
                assert time <= self.time
            else :
                print (f"Step {i} : No more events available !")
                print ("Stopping...")
                return
        #   Mets à jour l'efficacité quantique interne
        if self._recombinations > 0 :
            self._IQE = 100. * float(self._emission)/float(self._recombinations)
    


    ####################################
    ####____Méthodes get privées____####
    ####################################
    def _get_molecule_type(self, position : Point) -> type :
        return type(self._grid[position.z][position.y][position.x])

    def _get_molecule(self, position : Point) -> Host | TADF | Fluorescent :
        return self._grid[position.z][position.y][position.x]
    
    def _get_all_events(self) -> list[Event] :
        output : list[Event] = self._move_electron_events + self._move_hole_events \
        + self._move_exciton_events + self._isc_events + self._decay_events \
        + self._binding_events + self._capture_events
        return output
    


    ##################################################################
    ####____Méthode d'affichage du réseau dans son état actuel____####
    ##################################################################
    def plot(self, nom : str = "Inconnu") :
        use("Agg")
        plt.figure(dpi=100)
        axes = plt.axes(projection = "3d")

        x_grid = [0, self._dimension.x - 1]
        y_grid = [0, self._dimension.y - 1]
        z_grid = [0, self._dimension.z - 1]
        for x in x_grid :
            for y in y_grid :
                axes.plot([x,x], [y,y], z_grid, linestyle = "dashed", color = "k")
            for z in z_grid :
                axes.plot([x,x], y_grid, [z,z], linestyle = "solid", color = "k")
        for y in y_grid :
            for z in z_grid :
                axes.plot(x_grid, [y,y], [z,z], linestyle = "solid", color = "k")

        color = "b"
        electron_host = [
            position
            for position in self._electrons_locations
            if self._get_molecule_type(position) is Host
        ]
        if len(electron_host) > 0 :
            marker_style = "o"
            x = [position.x for position in electron_host]
            y = [position.y for position in electron_host]
            z = [position.z for position in electron_host]
            axes.scatter(x, y, z, s = 100, c = color, marker = marker_style)
        electron_tadf = [
            position
            for position in self._electrons_locations
            if self._get_molecule_type(position) is TADF
        ]
        if len(electron_tadf) > 0 :
            marker_style = "s"
            x = [position.x for position in electron_tadf]
            y = [position.y for position in electron_tadf]
            z = [position.z for position in electron_tadf]
            axes.scatter(x, y, z, s = 100, c = color, marker = marker_style)
        electron_fluorescent = [
            position
            for position in self._electrons_locations
            if self._get_molecule_type(position) is Fluorescent
        ]
        if len(electron_fluorescent) > 0 :
            marker_style = "^"
            x = [position.x for position in electron_fluorescent]
            y = [position.y for position in electron_fluorescent]
            z = [position.z for position in electron_fluorescent]
            axes.scatter(x, y, z, s = 75, c = color, marker = marker_style)

        color = "r"
        hole_host = [
            position
            for position in self._holes_locations
            if self._get_molecule_type(position) is Host
        ]
        if len(hole_host) > 0 :
            marker_style = "o"
            x = [position.x for position in hole_host]
            y = [position.y for position in hole_host]
            z = [position.z for position in hole_host]
            axes.scatter(x, y, z, s = 100, c = color, marker = marker_style)
        hole_tadf = [
            position
            for position in self._holes_locations
            if self._get_molecule_type(position) is TADF
        ]
        if len(hole_tadf) > 0 :
            marker_style = "s"
            x = [position.x for position in hole_tadf]
            y = [position.y for position in hole_tadf]
            z = [position.z for position in hole_tadf]
            axes.scatter(x, y, z, s = 100, c = color, marker = marker_style)
        hole_fluorescent = [
            position
            for position in self._holes_locations
            if self._get_molecule_type(position) is Fluorescent
        ]
        if len(hole_fluorescent) > 0 :
            marker_style = "^"
            x = [position.x for position in hole_fluorescent]
            y = [position.y for position in hole_fluorescent]
            z = [position.z for position in hole_fluorescent]
            axes.scatter(x, y, z, s = 100, c = color, marker = marker_style)

        color = "m"
        exciton_host = [
            position
            for position in self._excitons_locations
            if self._get_molecule_type(position) is Host
        ]
        if len(exciton_host) > 0 :
            marker_style = "o"
            x = [position.x for position in exciton_host]
            y = [position.y for position in exciton_host]
            z = [position.z for position in exciton_host]
            axes.scatter(x, y, z, s = 100, c = color, marker = marker_style)
        exciton_tadf = [
            position
            for position in self._excitons_locations
            if self._get_molecule_type(position) is TADF
        ]
        if len(exciton_tadf) > 0 :
            marker_style = "s"
            x = [position.x for position in exciton_tadf]
            y = [position.y for position in exciton_tadf]
            z = [position.z for position in exciton_tadf]
            axes.scatter(x, y, z, s = 100, c = color, marker = marker_style)
        exciton_fluorescent = [
            position
            for position in self._excitons_locations
            if self._get_molecule_type(position) is Fluorescent
        ]
        if len(exciton_fluorescent) > 0 :
            marker_style = "^"
            x = [position.x for position in exciton_fluorescent]
            y = [position.y for position in exciton_fluorescent]
            z = [position.z for position in exciton_fluorescent]
            axes.scatter(x, y, z, s = 100, c = color, marker = marker_style)

        axes.set_xlabel("x", size = 16)
        axes.set_ylabel("y", size = 16)
        axes.set_zlabel("z", size = 16)
        axes.set_xlim([0, self._dimension.x - 1])
        axes.set_ylim([0, self._dimension.y - 1])
        axes.set_zlim([0, self._dimension.z - 1])
        plt.tight_layout()
        plt.savefig(nom)
        plt.close()
        return

from time import time
start = time()
test = Lattice((10,10,10), (0.84, 0.15, 0.01), charges = 1)
op = 10**6
for i in range(100) :
    test.operations(1)
print(test._IQE)
end = time()
print(f"{end - start} s")