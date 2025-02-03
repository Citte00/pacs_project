#!/usr/bin/env python3

"""
@file errorvsize.py
@author Lorenzo Citterio (github.com/Citte00)
@date 2025-02-02

@copyright Copyright (c) 2025
"""

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy
import sys
import os

# Font.
matplotlib.rcParams.update({'font.size': 18})

# Sizes.
sizes: list[float] = []

# Degree (comparison).
degree: int = -1

# Errors.
l2_errors: list[float] = []
dg_errors: list[float] = []
energy_errors: list[float] = []

# File.
if len(sys.argv) <= 1:
    print(f"Usage: {sys.argv[0]} /path/to/file.")
    sys.exit(0)

try:
    if sys.argv[1].split(".")[-1] != "error":
        raise

    file = open(sys.argv[1], "r+")
    lines: list[str] = file.read().split("\n")
    file.close()

    # Title.
    title: str = lines[0]

except FileNotFoundError:
    print("File not found.")
    sys.exit(-1)

except:
    print("Load a .error file.")
    sys.exit(-1)

# Simple scraping.
for line in lines:
    try:
        data: list[str] = line.split(" ")

        if "Size" in line:
            sizes.append(float(data[-1]))

        elif "L2 Error" in line:
            l2_errors.append(float(data[-1]))

        elif "DG Error" in line:
            dg_errors.append(float(data[-1]))

        elif "Energy Error" in line:
            energy_errors.append(float(data[-1]))

        elif "Degree" in line:
            degree = int(data[-1])

    except ValueError:
        continue

# Colors.
black: list[float] = [7 / 255, 54 / 255, 66 / 255]
red: list[float] = [220 / 255, 50 / 255, 47 / 255]

# Comparison.
l2_comparison: list[float] = []
dg_comparison: list[float] = []
energy_comparison: list[float] = []

for size in sizes:
    l2_comparison.append((size / sizes[-1]) ** (degree + 1) * l2_errors[-1])
    dg_comparison.append((size / sizes[-1]) ** degree * dg_errors[-1])
    if energy_errors:
        energy_comparison.append((size / sizes[-1]) ** degree * energy_errors[-1])

# Ticks.
sizes_ticks = [sizes[0], sizes[-1]]

l2_ticks = [l2_errors[0], l2_errors[-1]]
dg_ticks = [dg_errors[0], dg_errors[-1]]
energy_ticks = [energy_errors[0], energy_errors[-1]] if energy_errors else []

# Labels.
sizes_labels = [f"{tick:.3f}" for tick in sizes_ticks]

l2_labels = [f"{tick:.1e}" for tick in l2_ticks]
dg_labels = [f"{tick:.1e}" for tick in dg_ticks]
energy_labels = [f"{tick:.1e}" for tick in energy_ticks] if energy_errors else []

# Plot.
fig, axes = plt.subplots(1, 3 if energy_errors else 2, figsize=(16,8))
if energy_errors:
    fig.suptitle("L2, DG and Energy errors vs. size")
else:
    fig.suptitle("L2 and DG errors vs. size")

# L2.
axes[0].plot(sizes, l2_errors, color=black, marker="*", linewidth=1, label="$L^2$ error") # Error.
axes[0].plot(sizes, l2_comparison, color=red, linestyle="--", linewidth=0.5, label=f"$h^{degree + 1}$") # Comparison.

# DG.
axes[1].plot(sizes, dg_errors, color=black, marker="*", linewidth=1, label="$DG$ error") # Error.
axes[1].plot(sizes, dg_comparison, color=red, linestyle="--", linewidth=0.5, label=f"$h^{degree}$") # Comparison.

# Energy error.
if energy_errors:
    axes[2].plot(sizes, energy_errors, color=black, marker="*", linewidth=1, label="Energy error")  # Error.
    axes[2].plot(sizes, energy_comparison, color=red, linestyle="--", linewidth=0.5, label=f"$h^{degree}$")  # Comparison.


if "lshape" in sys.argv[1]:
    dg_comparison: list[float] = []

    for size in sizes:
        dg_comparison.append((size / sizes[-1]) ** (2 / 3) * dg_errors[-1])

    axes[1].plot(sizes, dg_comparison, color=red, linestyle="--", linewidth=0.5, label="$h^{2/3}$") # Comparison.

# Parameters.
for j in range(2 + (1 if energy_errors else 0)):

    # Loglog scale.
    axes[j].set_xscale("log")
    axes[j].set_yscale("log")

    # Ticks.
    axes[j].xaxis.set_minor_formatter(NullFormatter())
    axes[j].yaxis.set_minor_formatter(NullFormatter())

    # Legend.
    axes[j].legend(loc="best")

axes[1].yaxis.tick_right()
axes[1].yaxis.set_label_position("right")

# Output.
if matplotlib.get_backend() == "agg":
    name = f"{os.path.splitext(sys.argv[1])[0]}.png"
    print(f"No interactive backend available. Plot will be saved as {name}.")
    plt.savefig(name)
else:
    plt.show()