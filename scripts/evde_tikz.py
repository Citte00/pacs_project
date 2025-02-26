#!/usr/bin/env python3

"""
@file evde_tikz.py
@brief errorvsize TikZ Python wrapper.
@author Lorenzo Citterio (github.com/Citte00)
@date 2025-02-26

@copyright Copyright (c) 2024
"""

import sys

# File.
if len(sys.argv) <= 1:
    print(f"Usage: {sys.argv[0]} /path/to/file.")
    sys.exit(0)

# Default template path
template_path = "templates/sizes_tikz.tex"

# Check if '--energy' argument is passed
if "--energy" in sys.argv:
    template_path = "templates/sizes_tikz_energy.tex"

# Template.
try:
    file = open(template_path, "r+")
    template: str = file.read()
    file.close()

except FileNotFoundError:
    print("Templates not found.")
    sys.argv(-1)

# Degree.
degrees: list[int] = []

# Errors.
l2_errors: list[float] = []
dg_errors: list[float] = []
energy_errors: list[float] = []

# Plots.
l2_plots: list[str] = []
dg_plots: list[str] = []
energy_plots: list[str] =[]

try:
    if sys.argv[1].split(".")[-1] != "error":
        raise

    file = open(sys.argv[1], "r+")
    lines: list[str] = file.read().split("\n")
    file.close()

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

        if "Degree" in line:
            degrees.append(int(data[-1]))

        elif "L2 Error" in line:
            l2_errors.append(float(data[-1]))

        elif "DG Error" in line:
            dg_errors.append(float(data[-1]))

        elif "Energy Error" in line:
            energy_errors.append(float(data[-1]))

    except ValueError:
        continue

# Plots.

# Errors.
l2_errors_string: str = "\n\\addplot[solarized-blue, mark=*] coordinates "
l2_errors_string += "{" + " ".join((f"({degrees[j]},{l2_errors[j]})" for j in range(len(degrees)))) + "};"
l2_errors_string += "\n"
l2_errors_string += "\\addlegendentry{$L^2$ Error}\n"

dg_errors_string: str = "\n\\addplot[solarized-blue, mark=*] coordinates "
dg_errors_string += "{" + " ".join((f"({degrees[j]},{dg_errors[j]})" for j in range(len(degrees)))) + "};"
dg_errors_string += "\n"
dg_errors_string += "\\addlegendentry{$DG$ Error}\n"

l2_plots.append(l2_errors_string)
dg_plots.append(dg_errors_string)

# Text.
template = template.replace("% PLOTS_L2", "".join(l2_plots))
template = template.replace("% PLOTS_DG", "".join(dg_plots))

# Energy error.
if "--energy" in sys.argv:
    
    energy_errors_string: str = "\n\\addplot[solarized-blue, mark=*] coordinates "
    energy_errors_string += "{" + " ".join((f"({degrees[j]},{energy_errors[j]})" for j in range(len(degrees)))) + "};"
    energy_errors_string += "\n"
    energy_errors_string += "\\addlegendentry{Energy Error}\n"
    energy_plots.append(energy_errors_string)

    template = template.replace("% PLOTS_ENERGY", "".join(energy_plots))

# Output.
file = open(sys.argv[1].replace(".error", "_evde.tex"), "w+")
file.write(template)
file.close()