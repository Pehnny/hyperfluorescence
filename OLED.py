#############################################################################################################
#
#   author(s) : Théo Piron
#   last update : 19/05/2022 (dd/mm/yyyy)
#   python version : 3.10.4
#   modules : reaseau, molecule, event, energy
#
#   commentaires :  Cette version exploite le multiprocessing. 
#                   Matplotlib n'est pas compatible avec le multithread.
#   
#
#############################################################################################################

from event import point
from reseau import lattice
from molecule import exciton, fluorescent, tadf, host
# from hote import host
# from TADF import tadf
# from emetteur import fluorescent

from datetime import datetime
from statistics import mean
from math import floor
import multiprocessing
from joblib import Parallel, delayed

###############################################################################################
#   Partie centrale du code : fait tourner les réseaux en parallèle
###############################################################################################

# Fonction qui permet de créer et d'utiliser les réseaux.
# Arguments :   - dimensions : dimensions du réseau (point)
#               - proportions : pourcentage de matériaux TADF et Fluo (tuple de taille 2)
#               - iterations : nombre d'étapes que l'on veut faire parcourir au réseau (int)
#               - moyenne : nombre de réseau que l'on veut faire tourner (int)
#               - plot : nombre d'étapes entre deux affichages du réseau (int) (default : None)
#                   - None => la fonction ne produira aucun plot
#               - Ez : intensité du champ électrique appliqué sur le réseau (float) (default : 10**8)
#               - sigma : désordre énergétique au sein des molécules du réseau (float) (default : 0.1)
def OLED(dimension : point | list[point], proportion : tuple | list[tuple], iterations : int, moyenne : int, plot : int = None, Ez : float = 10**8, sigma : float = 0.1) :
    date = str(datetime.today().strftime('%Y%m%d_%H%M%S'))
    size = str(dimension.x) + "x" + str(dimension.y) + "x" + str(dimension.z)
    name = date + "_" + size + "_" + str(iterations)
    # Réseau
    cores = floor(multiprocessing.cpu_count() / 2) # Nombre de coeurs utilisés 
    if __name__ == "__main__" :
        # Création des réseaux
        reseaux : list[lattice] = Parallel(n_jobs=cores, prefer="processes")(delayed(lambda dim, prop, E, hash, sig : lattice(dim, prop, E, hashtag=hash, sigma=sig))(dimension, proportion, Ez, i, sigma) for i in range(moyenne))
        # Proportions des types de réseaux
        # tadf_neighboors = []
        # for reseau in reseaux :
        #     for x in reseau.GRID :
        #         for y in x :
        #             for z in y :
        #                 if isinstance(z, tadf) :
        #                     fluo_neighboor = 0
        #                     for voisin in z.VOISINS :
        #                         if isinstance(reseau.GRID[voisin.x][voisin.y][voisin.z], fluorescent) :
        #                             fluo_neighboor += 1
        #                     tadf_neighboors.append(fluo_neighboor)
        # print("Nombre de voisins fluo aux tadf moyen : " + str(mean(tadf_neighboors)))
        # input("Waiting for any key")
        # Evolution du réseau
        if plot == None or plot > iterations :
            reseaux = Parallel(n_jobs=cores, prefer="processes")(delayed(_evolution)(reseau, iterations) for reseau in reseaux)
        elif plot < iterations :
            pauses = floor(iterations / plot)
            for i in range(pauses) :
                start = i * plot
                stop = (i+1) * plot
                for reseau in reseaux :
                    reseau.Plot(name + "_R" + str(reseau.HASHTAG) + "_" + str(start) + ".png")
                reseaux = Parallel(n_jobs=cores, prefer="processes")(delayed(_evolution)(reseau, stop, start, True) for reseau in reseaux)
            start = stop
            stop = iterations
            reseaux = Parallel(n_jobs=cores, prefer="processes")(delayed(_evolution)(reseau, stop, start, True) for reseau in reseaux)
            for reseau in reseaux :
                reseau.Plot(name + "_R" + str(reseau.HASHTAG) + "_" + str(stop) + ".png")        
        informations : list[infos] = Parallel(n_jobs=cores, prefer="processes")(delayed(_information)(reseau) for reseau in reseaux)
        # Positions des charges piégées
        _el = []
        _ho = []
        _el_host, _el_tadf, _el_fluo = [], [], []
        _ho_host, _ho_tadf, _ho_fluo = [], [], []
        for reseau in reseaux :
            for x in reseau.GRID :
                for y in x :
                    for z in y :
                        if z.electron :
                            _el.append(z)
                        if z.hole :
                            _ho.append(z)
            _hote, _tadf, _fluo = 0, 0, 0
            for m in _el :
                if isinstance(m, host) :
                    _hote += 1
                elif isinstance(m, tadf) :
                    _tadf += 1
                elif isinstance(m, fluorescent) :
                    _fluo += 1
            _el_host.append(_hote) 
            _el_tadf.append(_tadf)
            _el_fluo.append(_fluo)
            _hote, _tadf, _fluo = 0, 0, 0
            for m in _ho :
                if isinstance(m, host) :
                    _hote += 1
                elif isinstance(m, tadf) :
                    _tadf += 1
                elif isinstance(m, fluorescent) :
                    _fluo += 1
            _ho_host.append(_hote) 
            _ho_tadf.append(_tadf)
            _ho_fluo.append(_fluo)
            _el.clear()
            _ho.clear()
        print("electron sur hôte : " + str(mean(_el_host)))
        print("electron sur tadf : " + str(mean(_el_tadf)))
        print("electron sur fluo : " + str(mean(_el_fluo)))
        print("trou sur hôte : " + str(mean(_ho_host)))
        print("trou sur tadf : " + str(mean(_ho_tadf)))
        print("trou sur fluo : " + str(mean(_ho_fluo)))
    else :
        print ("OLED : processus non exécuté")
        return
    # Calcule les moyennes
    output = _moyenne(informations)
    # Informations générales 
    _ecriture(size, iterations, moyenne, output, name, proportion)
    # Ordre et localisation des excitons
    # name = name + "_excitons"
    # fichier = open(name + ".txt", "w")
    # count = 0
    # for i in reseau.where.Order :
    #     count += 1
    #     fichier.write(str(count) + " : " + i.__name__ + "\n")
    # fichier.close()
    print("Fin")
    return

###############################################################################################
#   Enregistre les informations intéressantes
###############################################################################################

# Classe qui stock les données pertinentes du réseau
class infos :
    def __init__(self, IQE : float, temps : float, electrons : float, trous : float, hostTot : float, tadfTot : float, fluoTot : float, E : float, sigma : float):
        self.IQE = IQE
        self.temps = temps
        self.el = electrons
        self.ho = trous
        self.Host = hostTot
        self.TADF = tadfTot
        self.Fluo = fluoTot
        self.E = E
        self.sigma = sigma

def _evolution(reseau : lattice, stop : int, start : int = 0, plot : bool = False) :
    if not plot :
        reseau.Operations(stop)
        return reseau
    else :
        reseau.Operations(stop, start)
        return reseau
# Fonction qui extrait les données pertinentes du réseau 
def _information(reseau : lattice) -> infos :
    IQE = reseau.IQE
    temps = reseau.temps
    electrons, trous = reseau.electrons, reseau.trous # excitons, électrons et trous présents dans le réseau à la fin de la simulation
    hostTot, tadfTot, fluoTot = reseau.where.HOST, reseau.where.TADF, reseau.where.FLUO # exctions totaux formés sur les hôtes, tadf et fluo
    E = reseau.E.z
    sigma = reseau.sigma
    return infos(IQE, temps, electrons, trous, hostTot, tadfTot, fluoTot, E, sigma)
# Fonction qui calcule la moyenne des données récoltées par infos
def _moyenne(data : infos | list[infos]) -> infos :
    IQE = []
    temps = []
    electrons, trous =  [], []
    hostTot, tadfTot, fluoTot = [], [], []
    E = []
    sigma = []
    if isinstance(data, list) :
        for i in data :
            if isinstance(i, infos) :
                IQE.append(i.IQE)
                temps.append(i.temps)
                electrons.append(i.el)
                trous.append(i.ho)
                hostTot.append(i.Host)
                tadfTot.append(i.TADF)
                fluoTot.append(i.Fluo)
                E.append(i.E)
                sigma.append(i.sigma)
            else :
                print ("moyenne : Erreur donnés de type inconnu")
                return
        return infos(mean(IQE), mean(temps), mean(electrons), mean(trous), mean(hostTot), mean(tadfTot), mean(fluoTot), mean(E), mean(sigma))
    elif isinstance(data, infos) :
        return data
# Fonction qui extrait les données dans un fichier
def _ecriture(size : str, iter : int, moy : int, info : infos, name : str, prop : tuple) -> None :
    fichier = open(name + ".txt", "w")
    string = "Reseau : " + size + "\n"
    string = string + "Etapes : " + str(iter) + "\n"
    string = string + "Moyenne : " + str(moy) + "\n"
    TADF, Fluo = prop
    Hote = 100 - TADF - Fluo
    string = string + "Proportions Host/TADF/Fluo : " + str(Hote) + "/" + str(TADF) + "/" + str(Fluo) + "\n"
    string = string + "IQE : " + str(info.IQE) + " %" + "\n"
    string = string + "Champ électrique : " + str(info.E) + " eV/m" + "\n"
    string = string + "Désordre énergétique : " + str(info.sigma) + " eV" + "\n"
    string = string + "Temps écoulé : " + str(info.temps) + " s" "\n"
    string = string + "Nombre d'électrons présents dans le réseau : " + str(info.el) + "\n"
    string = string + "Nombre de trous présents dans le réseau : " + str(info.ho) + "\n"
    string = string + "Nombre d'excitons formés sur les Host : " + str(info.Host) + "\n"
    string = string + "Nombre d'excitons formés sur les TADF : " + str(info.TADF) + "\n"
    string = string + "Nombre d'excitons formés sur les Fluo : " + str(info.Fluo) + "\n"
    # fichier.write("Nombre d'excitons présents sur les Host : " + str(host) + "\n")
    # fichier.write("Nombre d'excitons présents sur les TADF : " + str(tadf) + "\n")
    # fichier.write("Nombre d'excitons présents sur les Fluo : " + str(fluo) + "\n")
    # fichier.write(mol_type(reseau))
    fichier.write(string)
    fichier.close()
    return

###############################################################################################
#   Obsolète
###############################################################################################

def _mol_type(reseau : lattice) -> list :
    Host = 0
    Tadf = 0
    Fluo = 0
    for i in reseau.GRID :
        for j in i :
            for k in j :
                if reseau._isExciton(k.POS) :
                    if type(k) == host :
                        Host = Host + 1
                        continue
                    elif type(k) == tadf :
                        Tadf = Tadf + 1
                        continue
                    elif type(k) == fluorescent :
                        Fluo = Fluo + 1
                        continue
    return [Host, Tadf, Fluo]

def _how_many(reseau : lattice) :
    el = 0
    ho = 0
    ex = 0
    for i in reseau.GRID :
        for j in i :
            for k in j :
                if k.electron :
                    el += 1
                if k.hole :
                    ho += 1
                if isinstance(k.exciton, exciton) :
                    ex += 1
    if el != reseau.electrons :
        print("Erreur, pas le même nombre d'électrons")
    if ho != reseau.trous :
        print("Erreur, pas le même nombre d'électrons")
    if ex != reseau.excitons :
        print("Erreur, pas le même nombre d'électrons")

def where_charges(reseau : lattice) :
    output_el = []
    output_ho = []
    for i in reseau.GRID :
        for j in i :
            for k in j :
                if k.electron :
                    output_el.append(type(k))
                if k.hole :
                    output_ho.append(type(k))
                continue
    return output_el, output_ho

def print_where_charges(tab : list) :
    _hote = 0
    _tadf = 0
    _fluo = 0
    for i in tab :
        if isinstance(i, host) :
            _hote += 1
        elif isinstance(i, tadf) :
            _tadf += 1
        elif isinstance(i, fluorescent) :
            _fluo += 1
    return _hote, _tadf, _fluo

###############################################################################################
#   Exécution
###############################################################################################

OLED(point(10,10,5), (15, 1), 100, 1, plot = 1)
# # Référence
# OLED(point(15,15,10), (15, 1), 100000, 120, plot = 10000)
# Concentration des fluorescents
# OLED(point(15,15,10), (15, 5), 100000, 120, plot = 10000)
# # Concentration des TADF
# OLED(point(15,15,10), (5, 1), 100000, 100, plot = 10000)
# # Désordre énergétique
# OLED(point(15,15,10), (15, 1), 100000, 100, plot = 10000, sigma = 0)
# # Champ électrique
# OLED(point(15,15,10), (15, 1), 100000, 100, plot = 10000, Ez = 10**7)
# # Taille du réseau
# OLED(point(10,10,5), (15, 1), 100000, 100, plot = 10000)