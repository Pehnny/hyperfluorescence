from event import Point
from molecule import *
from matplotlib import pyplot as plt
from matplotlib import use

def plot(electrons : list[tuple[Point, type]], holes : list[tuple[Point, type]], excitons : list[tuple[Point, type]],
         x_size : int, y_size : int, z_size : int,
         name : str = "Inconnu") -> None :
    use("Agg")
    plt.figure(dpi=100)
    axes = plt.axes(projection = "3d")
    axes.set_xlabel("x", size = 16)
    axes.set_ylabel("y", size = 16)
    axes.set_zlabel("z", size = 16)
    axes.set_xlim([0, x_size - 1])
    axes.set_ylim([0, y_size - 1])
    axes.set_zlim([0, z_size - 1])

    x_grid = [0, x_size-1]
    y_grid = [0, y_size-1]
    z_grid = [0, z_size-1]
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
        for position, molecule in electrons
        if molecule is Host
    ]
    if len(electron_host) > 0 :
        marker_style = "o"
        x = [position.x for position in electron_host]
        y = [position.y for position in electron_host]
        z = [position.z for position in electron_host]
        axes.scatter(x, y, z, s = 100, c = color, marker = marker_style)
    electron_tadf = [
        position
        for position, molecule in electrons
        if molecule is TADF
    ]
    if len(electron_tadf) > 0 :
        marker_style = "s"
        x = [position.x for position in electron_tadf]
        y = [position.y for position in electron_tadf]
        z = [position.z for position in electron_tadf]
        axes.scatter(x, y, z, s = 100, c = color, marker = marker_style)
    electron_fluorescent = [
        position
        for position, molecule in electrons
        if molecule is Fluorophore
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
        for position, molecule in holes
        if molecule is Host
    ]
    if len(hole_host) > 0 :
        marker_style = "o"
        x = [position.x for position in hole_host]
        y = [position.y for position in hole_host]
        z = [position.z for position in hole_host]
        axes.scatter(x, y, z, s = 100, c = color, marker = marker_style)
    hole_tadf = [
        position
        for position, molecule in holes
        if molecule is TADF
    ]
    if len(hole_tadf) > 0 :
        marker_style = "s"
        x = [position.x for position in hole_tadf]
        y = [position.y for position in hole_tadf]
        z = [position.z for position in hole_tadf]
        axes.scatter(x, y, z, s = 100, c = color, marker = marker_style)
    hole_fluorescent = [
        position
        for position, molecule in holes
        if molecule is Fluorophore
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
        for position, molecule in excitons
        if molecule is Host
    ]
    if len(exciton_host) > 0 :
        marker_style = "o"
        x = [position.x for position in exciton_host]
        y = [position.y for position in exciton_host]
        z = [position.z for position in exciton_host]
        axes.scatter(x, y, z, s = 100, c = color, marker = marker_style)
    exciton_tadf = [
        position
        for position, molecule in excitons
        if molecule is TADF
    ]
    if len(exciton_tadf) > 0 :
        marker_style = "s"
        x = [position.x for position in exciton_tadf]
        y = [position.y for position in exciton_tadf]
        z = [position.z for position in exciton_tadf]
        axes.scatter(x, y, z, s = 100, c = color, marker = marker_style)
    exciton_fluorescent = [
        position
        for position, molecule in excitons
        if molecule is Fluorophore
    ]
    if len(exciton_fluorescent) > 0 :
        marker_style = "^"
        x = [position.x for position in exciton_fluorescent]
        y = [position.y for position in exciton_fluorescent]
        z = [position.z for position in exciton_fluorescent]
        axes.scatter(x, y, z, s = 100, c = color, marker = marker_style)

    plt.tight_layout()
    plt.savefig(name)
    plt.close()