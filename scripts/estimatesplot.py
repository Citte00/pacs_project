#!/usr/bin/env python3

"""
@file batch_estimatesplot.py
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
from glob import glob

# Font.
matplotlib.rcParams.update({'font.size': 18})

def process_poly_file(file_path):
    """Process a single .poly file and generate a .png file."""
    try:
        with open(file_path, "r") as file:
            lines = file.read().split("\n")

    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return

    except:
        print(f"Error reading file: {file_path}")
        return

    # Black.
    black = [7 / 255, 54 / 255, 66 / 255]

    # Normalize the colormap
    estimates = []
    errors = []
    comparison = []

    for line in lines:
        if line and line[0] != "@":
            data = line.split(" ")
            try:
                estimates.append(float(data[-3]))
                errors.append(float(data[-2]))
                comparison.append(float(data[-1]))
            except ValueError:
                continue

    if not (estimates and errors and comparison):
        print(f"No valid data in {file_path}")
        return
    
    estimates_list = [estimates, errors, comparison]

    # Define titles for each plot
    titles = ["Estimates", "DG error", "Comparison"]

    # Create a figure and axis
    fig, ax = plt.subplots(1, 3, figsize=(24, 8))
    fig.suptitle("Error Estimates Distribution", color="white", fontsize=20)

    for k in range(len(estimates_list)):

        norm = Normalize(vmin=min(estimates_list[k]), vmax=max(estimates_list[k]))

        for line in lines:
            if line and line[0] != "@":
                x, y = [], []
                data = [float(number) for number in line.split(" ") if number]

                try:
                    for j in range(0, len(data) - 3, 2):
                        x.append(data[j])
                        y.append(data[j + 1])

                except ValueError:
                    continue

                if not (x and y):
                    continue

                # Color.
                color = list(cm.turbo(norm(float(data[-3 + k]))))
                color[3] = 0.75  # Reduces alpha.

                # Plot.
                ax[k].fill(x, y, facecolor=color, edgecolor=black, linewidth=0.25)

        # Create a ScalarMappable for the colorbar
        sm = plt.cm.ScalarMappable(cmap=cm.turbo, norm=norm)
        sm.set_array([])

        # Add colorbar
        cbar = fig.colorbar(sm, ax=ax[k], ticks=np.linspace(min(estimates_list[k]), max(estimates_list[k]), num=5), alpha=0.75)
        cbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))

        ax[k].set_title(titles[k], color="white", fontsize=18)
        ax[k].set_aspect('equal', adjustable='box')

        # Background and axes styling
        fig.patch.set_facecolor("#3A3A3A")
        ax[k].set_facecolor("#3A3A3A")

        # Set text and borders to light gray/white for contrast
        for spine in ax[k].spines.values():
            spine.set_color("#D3D3D3")

        ax[k].xaxis.label.set_color("white")
        ax[k].yaxis.label.set_color("white")
        ax[k].tick_params(colors="white")

        cbar.ax.yaxis.label.set_color("white")
        cbar.ax.tick_params(colors="white")

    # Save the output in the same directory as the input file
    output_name = f"{os.path.splitext(file_path)[0]}.png"
    plt.savefig(output_name)
    print(f"Saved: {output_name}")
    plt.close(fig)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} /path/to/file_or_folder")
        sys.exit(0)

    input_path = sys.argv[1]

    if os.path.isfile(input_path) and input_path.endswith(".poly"):
        # Process a single file
        process_poly_file(input_path)

    elif os.path.isdir(input_path):
        # Process all .poly files in the directory
        poly_files = glob(os.path.join(input_path, "*.poly"))

        if not poly_files:
            print(f"No .poly files found in {input_path}.")
            sys.exit(0)

        for poly_file in poly_files:
            process_poly_file(poly_file)

    else:
        print(f"Error: {input_path} is not a valid file or directory.")
        sys.exit(1)