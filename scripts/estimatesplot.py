#!/usr/bin/env python3

"""
@file estimatesplot.py
@author Lorenzo Citterio (github.com/Citte00)
@date 2025-01-24

@copyright Copyright (c) 2024
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import cm
from matplotlib.colors import Normalize
import numpy as np
import sys
import os

# Font.
matplotlib.rcParams.update({'font.size': 18})

if len(sys.argv) <= 1:
    print(f"Usage: {sys.argv[0]} /path/to/file.poly.")
    sys.exit(0)

try:
    if sys.argv[1].split(".")[-1] != "poly":
        raise

    file = open(sys.argv[1], "r+")
    lines: list[str] = file.read().split("\n")
    file.close()

except FileNotFoundError:
    print("File not found.")
    sys.exit(-1)

except:
    print("Load a .poly file.")
    sys.exit(-1)

# Black.
black: list[float] = [7 / 255, 54 / 255, 66 / 255]

# Create a figure and axis
fig, ax = plt.subplots(figsize=(10,8))

# Normalize the colormap
estimates: list[float] = []


for line in lines:
    if line:
        if line[0] == "@":
            continue

    data: list[str] = line.split(" ")
    
    try:
        estimates.append(float(data[-1]))

    except ValueError:
        continue

norm = Normalize(vmin=min(estimates), vmax=max(estimates))

for line in lines:
    if line:
        if line[0] == "@":
            continue

    x: list[float] = []
    y: list[float] = []

    data: list[str] = line.split(" ")
    data: list[float] = [float(number) for number in data if number]
    
    try:
        for j in range(0, len(data) if len(data) % 2 == 0 else len(data) - 1, 2):
            x.append(float(data[j]))
            y.append(float(data[j + 1]))

    except ValueError:
        continue

    if not (x and y):
        continue

    # Color.
    color: tuple[int] = list(cm.Blues(norm(float(data[-1]))))
    color[3] = 0.75 # Reduces alpha.

    # Plot.
    ax.fill(x, y, facecolor=color, edgecolor=black, linewidth=0.25)

# Create a ScalarMappable for the colorbar
sm = plt.cm.ScalarMappable(cmap=cm.Blues, norm=norm)
sm.set_array([])

# Add colorbar
cbar = fig.colorbar(sm, ax=ax, ticks=np.linspace(min(estimates), max(estimates), num=5), alpha=0.75)
cbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))  # Adjust precision as needed

ax.set_aspect('equal', adjustable='box')

# Output.
if matplotlib.get_backend() == "agg":
    name = f"{os.path.splitext(sys.argv[1])[0]}.png"
    print(f"No interactive backend available. Plot will be saved as {name}.")
    plt.savefig(name)
else:
    plt.show()