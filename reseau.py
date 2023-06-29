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
from math import exp, log, prod, inf
from random import Random
from collections import deque

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
#         self.Order.append(Fluorophore)

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
    """Classe représentant un réseau cristallin de type OLED hyperFluorophoree.

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
    _charge_transfer_rate : float
        Taux de transfert des charges au sein du réseau.
    _temperature : float
        Température de fonctionnement du réseau.
    _grid : list[list[list[Host | TADF | Fluorophore]]]
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
    _lattice_creation() -> list[list[list[Host | TADF | Fluorophore]]]
        ...
    _molecule_type(n : int, position : Point) -> Host | TADF | Fluorophore
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
    def __init__(self, proportions : tuple[float,float,float], dimension : tuple[int,int,int] = (20, 20, 10),
                 electric_field : float = 10.**(-1), charges : int = 4, charge_tranfer_distance : int = 1,
                 cutoff_radius : float = 19.2, architecture : str = NotImplemented) -> None :
        self._init_raises(dimension, proportions, charge_tranfer_distance)
        self._seed : Random = Random()
        self._lattice_parameters_creation(dimension, proportions, electric_field, charges, cutoff_radius)
        self._grid : list[list[list[Host | TADF | Fluorophore]]] = self._lattice_creation(charge_tranfer_distance)
        self._charges_injection()
        self._events_creation()
        self._injection : int = 2 * charges
        self._emission : int = 0
        self._recombination : int = 0
        self._IQE : float = 0.
        self._step : int = 0
        self._time : float = 0.
        self._cache : deque[Event] = deque((None for i in range(10)), 10)

    def _init_raises(self, dimension : tuple[int, int, int], proportions : tuple[float, float, float],
                     distance : int) -> None :
        if distance >= min(dimension) :
            raise ValueError(f"Expected charge_transfer_distance to be lower than all 3 dimensions, got {distance}.")
        if not isinstance(dimension, tuple) :
            raise TypeError(f"Expected type(dimension) to be tuple, got {type(dimension)}.")
        elif len(dimension) != 3 :
            raise IndexError(f"Expected dimension length to be 3, got {len(dimension)}.")
        elif dimension[-1] < 3 :
            raise ValueError(f"Expected the last component of dimension to be > 2, got {dimension[-1]}.")
        if not isinstance(dimension, tuple) :
            raise TypeError(f"Expected type(proportions) to be tuple, got {type(dimension)}.")
        elif len(proportions) != 3 :
            raise IndexError(f"Expected proportions length to be 3, got {len(dimension)}.")
        if  min(proportions) * prod(dimension) < 1 :
            minimum = 1 / min(proportions)
            raise ValueError(f"Not enough molecules for current proportions and dimension.\n Expected at least {minimum}, got {prod(dimension)}.")
        
    def _lattice_parameters_creation(self, dimension : tuple[int,int,int], proportions : tuple[float,float,float],
                                     electric_field : float, charges : int, cutoff_radius : float) -> None :
        if sum(proportions) != 1. :
            norm = sum(dimension)
            proportions = (dimension[0] / norm, dimension[1] / norm, dimension[2] / norm)
        self._dimension : Point = Point(*dimension)
        self._proportions : Proportion = Proportion(*proportions)
        self._electric_field : Vector = Vector(0, 0, electric_field)    # [eV/nm]
        self._lattice_constant : float = 1.                             # [nm]
        self._gamma : float = 10. / self._lattice_constant              # [1/nm]
        self._temperature : float = 300.                                # [K]
        self._cutoff : float = cutoff_radius                            # [nm]
        self._charges : int = charges

    def _lattice_creation(self, distance) -> list[list[list[Host | TADF | Fluorophore]]] :
        self._fluorophores_locations : list[Point] = []
        x_max : int = self._dimension.x
        y_max : int = self._dimension.y
        z_max : int = self._dimension.z
        grid_size : int = x_max * y_max * z_max
        host_layers : int = int(grid_size * self._proportions.host) // (x_max * y_max)
        if host_layers%2 :
            host_layers -= 1
        elif host_layers > 0 :
            host_layers -= 2
        n_tadf : int = round(grid_size * self._proportions.tadf)
        n_fluo : int = round(grid_size * self._proportions.fluo)
        n_host : int = grid_size - host_layers * x_max * y_max - n_tadf - n_fluo
        assert n_host+n_tadf+n_fluo+host_layers*x_max*y_max == grid_size, f"Expected {grid_size} molecules, got {n_fluo+n_host+n_tadf}."
        sub_z_max : int = z_max - host_layers
        sub_grid : list[int] = self._seed.sample(
            [0, 1, 2],
            k = n_host + n_tadf + n_fluo,
            counts = [n_host, n_tadf, n_fluo]
        )
        assert len(sub_grid) == x_max * y_max * sub_z_max, f"Expected len(sub_grid) to be {x_max * y_max * sub_z_max}, got {len(sub_grid)}."
        grid : list[list[list[int]]] = [[[0 for x in range(x_max)] for y in range(y_max)] for z in range(host_layers//2)]
        grid.extend([[sub_grid[y * x_max : (y+1) * x_max] for y in range(y_max)] for z in range(sub_z_max)])
        grid.extend([[[0 for x in range(x_max)] for y in range(y_max)] for z in range(host_layers//2)])
        return [[[self._molecule_type(n, Point(x,y,z), distance) for x, n in enumerate(ssgrid)] for y, ssgrid in enumerate(sgrid)] for z, sgrid in enumerate(grid)]
    
    def _molecule_type(self, n : int, position : Point, distance : int) -> Host | TADF | Fluorophore :
        if n == 0 :
            return Host(position, self._neighbourhood(position, distance))
        elif n == 1 :
            return TADF(position, self._neighbourhood(position, distance))
        elif n == 2 :
            self._fluorophores_locations.append(position)
            return Fluorophore(position, self._neighbourhood(position, distance))
        raise ValueError(f"n should be 0, 1 or 2, got {n}")

    def _neighbourhood(self, position : Point, distance : int) -> list[Point] :
        x_range = self._born_von_karman(position.x, distance, "x")
        y_range = self._born_von_karman(position.y, distance, "y")
        z_range = self._not_born_von_karman(position.z, distance)
        return [Point(x,y,z) for x in x_range for y in y_range for z in z_range if (x, y, z) != (position.x, position.y, position.z)]

    def _born_von_karman(self, position : int, distance : int, axe : str) -> list[int] :
        size : int = getattr(self._dimension, axe)
        lower_bound : int = position - distance
        upper_bound : int = position + distance
        if lower_bound > -1 and upper_bound < size :
            return list(range(position - distance, position + distance + 1))
        elif lower_bound < 0 :
            return list(range(distance + 1)) + list(range(size - distance, size))
        else :
            return list(range(position - distance, size)) + list(range(distance))
    
    def _not_born_von_karman(self, position : int, distance : int) -> list[int] :
        size : int = self._dimension.z
        lower_bound : int = position - distance
        upper_bound : int = position + distance
        if lower_bound > -1 and upper_bound < size :
            return list(range(position - distance, position + distance + 1))
        elif lower_bound < 0 :
            return list(range(distance + 1))
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
        self._electron_events : list[Event] = self._init_move_electron_events()
        self._hole_events : list[Event] = self._init_move_hole_events()
        self._exciton_events : list[Event] = []

    def _init_move_electron_events(self) -> list[Event] :
        molecules = (
            self._get_molecule(position)
            for position in self._electrons_locations
        )
        events = (
            [
                Event(molecule.position, neighbour, self._time_move_electron(molecule.position, neighbour), EVENTS["move"], PARTICULES["electron"])
                for neighbour in molecule.neighbourhood
                if neighbour not in self._electrons_locations
            ]
            for molecule in molecules
        )
        fastests = [
            min(event)
            for event in events
        ]
        return fastests

    def _init_move_hole_events(self) -> list[Event] :
        molecules = (
            self._get_molecule(position)
            for position in self._holes_locations
        )
        events = (
            [
                Event(molecule.position, neighbour, self._time_move_hole(molecule.position, neighbour), EVENTS["move"], PARTICULES["hole"])
                for neighbour in molecule.neighbourhood
                if neighbour not in self._holes_locations
            ]
            for molecule in molecules
        )
        fastests = [
            min(event)
            for event in events
        ]
        return fastests



    ##################################################################
    ####____Méthodes de calcul des durées de chaque événement_____####
    ##################################################################
    def _time_move_electron(self, initial : Point, final : Point) -> float :
        rng : float = 1. - self._seed.random()
        movement : Vector = (final - initial) * self._lattice_constant
        delta_energy : float = self._lumo_energy(initial, final)
        delta_energy += 1. * self._electric_field * movement
        delta_energy += self._electron_electrostatic_energy(initial, final)
        transfer_rate : float = TRANSFER_RATES["charges"] * exp(-2. * self._gamma * movement.norm())
        if delta_energy >= 0 :
            transfer_rate *= exp(- delta_energy / (cst.BOLTZMANN * self._temperature))
        return - log(rng) / transfer_rate
        
    def _lumo_energy(self, initial : Point, final : Point) -> float :
        return self._get_molecule(final).lumo_energy - self._get_molecule(initial).lumo_energy
    
    def _electron_electrostatic_energy(self, initial : Point, final : Point) -> float :
        output : float = 0.
        for location in self._holes_locations :
            if location == final :
                return -inf
            old_ir : float = self._inverse_radius(((location - initial) * self._lattice_constant).norm())
            new_ir : float = self._inverse_radius(((location - final) * self._lattice_constant).norm())
            delta_ir : float = new_ir - old_ir
            output -= cst.ELECTROSTATIC * delta_ir
        for location in self._electrons_locations :
            if location == initial : continue
            old_ir : float = self._inverse_radius(((location - initial) * self._lattice_constant).norm())
            new_ir : float = self._inverse_radius(((location - final) * self._lattice_constant).norm())
            delta_ir : float = new_ir - old_ir
            output += cst.ELECTROSTATIC * delta_ir
        return output

    def _time_move_hole(self, initial : Point, final : Point) -> float :
        rng = 1. - self._seed.random()
        movement : Vector = (final - initial) * self._lattice_constant
        delta_energy = self._homo_energy(initial, final)
        delta_energy += -1. * self._electric_field * movement
        delta_energy += self._hole_electrostatic_energy(initial, final)
        transfer_rate : float = TRANSFER_RATES["charges"] * exp(-2. * self._gamma * movement.norm())
        if delta_energy >= 0 :
            transfer_rate *= exp(- delta_energy / (cst.BOLTZMANN * self._temperature))
        return - log(rng) / transfer_rate
        
    def _homo_energy(self, initial : Point, final : Point) -> float :
        return self._get_molecule(final).homo_energy - self._get_molecule(initial).homo_energy

    def _hole_electrostatic_energy(self, initial : Point, final : Point) -> float :
        output : float = 0.
        for location in self._electrons_locations :
            if location == final :
                return -inf
            old_ir : float = self._inverse_radius(((location - initial) * self._lattice_constant).norm())
            new_ir : float = self._inverse_radius(((location - final) * self._lattice_constant).norm())
            delta_ir : float = new_ir - old_ir
            output -= cst.ELECTROSTATIC * delta_ir
        for location in self._holes_locations :
            if location == initial : continue
            old_ir : float = self._inverse_radius(((location - initial) * self._lattice_constant).norm())
            new_ir : float = self._inverse_radius(((location - final) * self._lattice_constant).norm())
            delta_ir : float = new_ir - old_ir
            output += cst.ELECTROSTATIC * delta_ir
        return output
    
    def _inverse_radius(self, distance : float) -> float :
        if distance > self._cutoff :
            return 0.
        else : 
            return 1. / distance - 1. / self._cutoff
    
    def _time_FRET(self, initial : Point, final : Point,
                   spin : int) -> float :
        distance : float = ((final - initial) * self._lattice_constant).norm()
        rng : float = 1. - self._seed.random()
        if spin == EXCITON["triplet"] :
            transfer_rate : float = (TRANSFER_RADIUS["TTS"] / distance)**6 * TRANSFER_RATES["ACRSA_F"]
        elif spin == EXCITON["singlet"] :
            transfer_rate : float = (TRANSFER_RADIUS["STS"] / distance)**6 * TRANSFER_RATES["ACRSA_F"]
        return - log(rng) / transfer_rate
        
    def _time_ISC(self, position : Point, spin : int) -> float :
        rng : float = 1. - self._seed.random()
        delta_EST : float = self._singlet_to_triplet_energy(position)
        if spin == EXCITON["triplet"] :
            transfer_rate : float = SOC["RISC"] * exp(- delta_EST / (cst.BOLTZMANN * self._temperature))
        elif spin == EXCITON["singlet"] :
            transfer_rate : float = SOC["ISC"] * exp(delta_EST / (cst.BOLTZMANN * self._temperature))
        return - log(rng) / transfer_rate
    
    def _singlet_to_triplet_energy(self, position : Point) -> float :
        molecule = self._get_molecule(position)
        return molecule.s1_energy - molecule.t1_energy
    
    def _time_TADF_decay(self, spin : int) -> float :
        rng : float = 1. - self._seed.random()
        if spin == EXCITON["triplet"] :
            transfer_rate : float = TRANSFER_RATES["ACRSA_PH"]
        elif spin == EXCITON["singlet"] :
            transfer_rate : float = TRANSFER_RATES["ACRSA_F"]
        return - log(rng) / transfer_rate
        
    def _time_fluorophore_decay(self, spin : int) -> float :
        rng : float = 1. - self._seed.random()
        if spin == EXCITON["triplet"] :
            transfer_rate : float = TRANSFER_RATES["TBPe_NR"]
        elif spin == EXCITON["singlet"] :
            transfer_rate : float = TRANSFER_RATES["TBPe_F"]
        return - log(rng) / transfer_rate


    
    ################################################################################
    ####____Méthodes qui suppriment les événements qui ne sont plus utilisés____####
    ################################################################################
    def _remove_move_electron_events(self, event : Event) -> None :
        removed_events : list[Event] = []
        while True :
            try :
                index = self._electron_events.index(event)
                removed_events.append(self._electron_events[index])
                self._electron_events.remove(event)
            except ValueError :
                if len(removed_events) > 1 :
                    for item in removed_events :
                        if item.initial != event.initial :
                            self._new_move_electron_event(item.initial)
                return
        
    def _remove_move_hole_events(self, event : Event) -> None :
        removed_events : list[Event] = []
        while True :
            try :
                index = self._hole_events.index(event)
                removed_events.append(self._hole_events[index])
                self._hole_events.remove(event)
            except ValueError :
                if len(removed_events) > 1 :
                    for item in removed_events :
                        if item.initial != event.initial :
                            self._new_move_hole_event(item.initial)
                break

    def _remove_bound_event(self, event : Event) -> None :
        self._exciton_events.remove(event)

    def _remove_decay_event(self, event : Event) -> None :
        self._exciton_events.remove(event)

    def _remove_capture_electron_event(self, event : Event) -> None :
        self._electron_events.remove(event)

    def _remove_capture_hole_event(self, event : Event) -> None :
        self._hole_events.remove(event)

    def _remove_ISC_event(self, event : Event) -> None :
        self._exciton_events.remove(event)



    ############################################################################
    ####____Méthodes qui génèrent les nouveaux événements à chaque étape____####
    ############################################################################
    def _new_move_electron_event(self, position : Point) -> None :
        neighbourhood = self._get_molecule(position).neighbourhood
        events = [
            Event(position, neighbour, self._time_move_electron(position, neighbour), EVENTS["move"], PARTICULES["electron"])
            for neighbour in neighbourhood
            if not self._get_molecule(neighbour).electron
        ]
        self._electron_events.append(min(events))

    def _new_move_hole_event(self, position : Point) -> None :
        neighbourhood = self._get_molecule(position).neighbourhood
        events = [
            Event(position, neighbour, self._time_move_hole(position, neighbour), EVENTS["move"], PARTICULES["hole"])
            for neighbour in neighbourhood
            if not self._get_molecule(neighbour).hole
        ]
        self._hole_events.append(min(events))

    def _new_bound_event(self, position : Point) -> None :
        event = Event(position, position, 0., EVENTS["bound"], PARTICULES["exciton"])
        self._exciton_events.append(event)

    def _new_capture_electron_event(self, position : Point) -> None :
        event = Event(position, position, 0., EVENTS["capture"], PARTICULES["electron"])
        self._electron_events.append(event)

    def _new_capture_hole_event(self, position : Point) -> None :
        event = Event(position, position, 0., EVENTS["capture"], PARTICULES["hole"])
        self._hole_events.append(event)

    def _new_host_decay_event(self, position : Point) -> None :
        event = Event(position, position, 0., EVENTS["decay"], PARTICULES["exciton"])
        self._exciton_events.append(event)

    def _new_TADF_ISC_event(self, position : Point, spin : int) -> Event :
        event = Event(position, position, self._time_ISC(position, spin), EVENTS["ISC"], PARTICULES["exciton"])
        return event
    
    def _new_TADF_decay_event(self, position : Point, spin : int) -> Event :
        event = Event(position, position, self._time_TADF_decay(spin), EVENTS["decay"], PARTICULES["exciton"])
        return event

    def _new_TADF_FRET_event(self, position : Point, spin : int) -> Event :
        events = [
            Event(position, transfer, self._time_FRET(position, transfer, spin), EVENTS["Forster"], PARTICULES["exciton"])
            for transfer in self._fluorophores_locations
        ]
        event : Event = min(list(events))
        return event
    
    def _new_TADF_event(self, event : Event) -> None :
        self._exciton_events.append(event)

    def _new_fluorophore_decay_event(self, position : Point, spin : int) -> None :
        event = Event(position, position, self._time_fluorophore_decay(spin), EVENTS["decay"], PARTICULES["exciton"])
        self._exciton_events.append(event)



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
        self._injection += 1

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
        self._injection += 1

    def _decay(self, position : Point) -> None :
        photon = self._grid[position.z][position.y][position.x].exciton_decay()
        self._excitons_locations.remove(position)
        self._recombination += 1
        if photon : self._emission += 1

    def _intersystem_crossing(self, position : Point) -> None :
        self._grid[position.z][position.y][position.x].intersystem_crossing()

    def _FRET(self, initial : Point, final : Point) -> None :
        self._grid[initial.z][initial.y][initial.x].exciton_decay()
        self._excitons_locations.remove(initial)
        molecule = self._grid[final.z][final.y][final.x]
        molecule.exciton = EXCITON["singlet"]
        photon = molecule.exciton_decay()
        self._recombination += 1
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
        self._cache.popleft()
        self._cache.append(event)
        #   Traite les événements de type "move"
        if event.kind == EVENTS["move"] :
            #   Traite le cas d'un électron
            if event.particule == PARTICULES["electron"] :
                self._move_electron(event.initial, event.final)
                self._remove_move_electron_events(event)
                molecule = self._get_molecule(event.final)
                if not molecule.hole and event.final.z != 0 :
                    self._update_events(event.tau)
                    self._new_move_electron_event(event.final)
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
                self._move_hole(event.initial, event.final)
                self._remove_move_hole_events(event)
                molecule = self._get_molecule(event.final)
                if not molecule.electron and event.final.z != (self._dimension.z - 1) :
                    self._update_events(event.tau)
                    self._new_move_hole_event(event.final)
                elif molecule.electron :
                    event = Event(event.final, event.final, 0., EVENTS["move"], PARTICULES["electron"])
                    self._remove_move_electron_events(event)
                    self._update_events(event.tau)
                    self._new_bound_event(event.final)
                elif event.final.z == (self._dimension.z - 1) :
                    self._update_events(event.tau)
                    self._new_capture_hole_event(event.final)
        #   Traite les événements de type formation d'exciton
        elif event.kind == EVENTS["bound"] :
            self._form_exciton(event.final)
            self._remove_bound_event(event)
            molecule = self._get_molecule(event.final)
            if isinstance(molecule, Host) :
                self._new_host_decay_event(event.final)
            elif isinstance(molecule, TADF) :
                TADF_events : list[Event] = [
                    self._new_TADF_decay_event(event.final, molecule.exciton),
                    self._new_TADF_FRET_event(event.final, molecule.exciton),
                    self._new_TADF_ISC_event(event.final, molecule.exciton)
                ]
                self._new_TADF_event(min(TADF_events))
            elif isinstance(molecule, Fluorophore) :
                self._new_fluorophore_decay_event(event.final, molecule.exciton)
        #   Traite les événements de type conversion intersystème
        elif event.kind == EVENTS["ISC"] :
            self._intersystem_crossing(event.initial)
            self._remove_ISC_event(event)
            self._update_events(event.tau)
            molecule = self._get_molecule(event.final)
            decay_events : list[Event] = [
                    self._new_TADF_decay_event(event.final, molecule.exciton),
                    self._new_TADF_FRET_event(event.final, molecule.exciton),
            ]
            self._new_TADF_event(min(decay_events))
        #   Traite les événements de type transfert d'énergie de Forster
        elif event.kind == EVENTS["Forster"] :
            self._decay(event.initial)
            self._remove_decay_event(event)
            self._update_events(event.tau)
            self._electron_reinjection()
            self._new_move_electron_event(self._electrons_locations[-1])
            self._hole_reinjection()
            self._new_move_hole_event(self._holes_locations[-1])
        #   Traite les événements de type recombinaison
        elif event.kind == EVENTS["decay"] :
            self._decay(event.initial)
            self._remove_decay_event(event)
            self._update_events(event.tau)
            self._electron_reinjection()
            self._new_move_electron_event(self._electrons_locations[-1])
            self._hole_reinjection()
            self._new_move_hole_event(self._holes_locations[-1])
        #   Traite les événements de type séparation d'exciton :
        elif event.kind == EVENTS["capture"] :
            if event.particule == PARTICULES["electron"] :
                self._capture_electron(event.final)
                self._remove_capture_electron_event(event)
                self._electron_reinjection()
                self._new_move_electron_event(self._electrons_locations[-1])
            elif event.particule == PARTICULES["hole"] :
                self._capture_hole(event.final)
                self._remove_capture_hole_event(event)
                self._hole_reinjection()
                self._new_move_hole_event(self._holes_locations[-1])
        return True
    
    def _update_events(self, time : float) -> None :
        for event in self._electron_events :
            event.tau -= time
        for event in self._hole_events :
            event.tau -= time
        for event in self._exciton_events :
            event.tau -= time

    def operations(self, recombinations : int, stop : int = 10**8) -> None :
        count : int = 0
        while self._recombination < recombinations :
            count += 1
            self._step += 1
            time = self._time
            #   Exécute l'évenement suivant et s'assure que le temps n'a pas diminué.
            try : 
                running = self._first_reaction_method()
            except ZeroDivisionError :
                print("Zero Division occured.\n", "Stopping process...") 
                print(self._cache)
                return
            if not running :
                print("No more events left.\n", "Stopping process...")
                self._IQE = 100. * 2. * float(self._emission) / float(self._injection)
                return 
            if time > self._time :
                print("Negative time occured.\n", "Stopping process...")
                return
            if count == stop :
                # print("Occurrence limit reached.\n", "Stopping process...")
                self._IQE = 100. * 2. * float(self._emission) / float(self._injection)
                return
        print("Recquired amount of recombinations reached.\n", "Stopping process...") 
        self._IQE = 100. * 2. * float(self._emission) / float(self._injection)
    


    ####################################
    ####____Méthodes get____####
    ####################################
    def _get_molecule_type(self, position : Point) -> type :
        return type(self._grid[position.z][position.y][position.x])

    def _get_molecule(self, position : Point) -> Host | TADF | Fluorophore :
        return self._grid[position.z][position.y][position.x]
    
    def _get_all_events(self) -> list[Event] :
        output : list[Event] = self._electron_events + self._hole_events + self._exciton_events
        return output
    
    def get_IQE(self) -> float :
        return self._IQE
    
    def get_particules_positions(self) -> tuple[list[tuple[Point, type]], list[tuple[Point, type]], list[tuple[Point, type]]] :
        electrons_locations = [(position, self._get_molecule_type(position)) for position in self._electrons_locations]
        holes_locations = [(position, self._get_molecule_type(position)) for position in self._holes_locations]
        excitons_locations = [(position, self._get_molecule_type(position)) for position in self._excitons_locations]
        return electrons_locations, holes_locations, excitons_locations
    
    def get_dimensions(self) -> tuple[int, int, int] :
        return self._dimension.x, self._dimension.y, self._dimension.z
