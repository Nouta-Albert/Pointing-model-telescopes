# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 10:48:40 2025

@author: Albert Einstein
"""

import numpy as np
import matplotlib.pyplot as plt

def spiral_grid(n, ntour=5):
    """Generate a simple Archimedean spiral (r, theta)."""
    th = np.linspace(0, ntour * 2 * np.pi, n)
    r = th / th.max()
    return r, th

def biased_fibonacci_grid(n, bias=0.8):
    """Generate biased Fibonacci spiral toward zenith using golden angle."""
    golden_angle = np.pi * (3 - np.sqrt(5))
    idx = np.arange(n)
    r = np.linspace(0, 1, n) ** bias  # bias < 1 → denser near zenith
    th = idx * golden_angle
    return r, th

def map_to_alt_az(r, th, alt_min_deg=30, alt_max_deg=85):
    """Map normalized radius to true altitude and azimuth."""
    alt = alt_max_deg - r * (alt_max_deg - alt_min_deg)
    az = np.degrees(th) % 360
    return alt, az

def plot_grids_with_altitude(grid_sizes=[30, 60, 100], alt_min=30, alt_max=85, bias=0.8):
    fig, axs = plt.subplots(len(grid_sizes), 2, subplot_kw={'projection': 'polar'}, figsize=(12, 10))
    fig.suptitle(f"Spiral vs Biased Fibonacci (bias = {bias}) Sampling", fontsize=16)

    for i, n in enumerate(grid_sizes):
        # Spiral
        r1, th1 = spiral_grid(n)
        alt1, az1 = map_to_alt_az(r1, th1, alt_min, alt_max)
        axs[i, 0].scatter(np.radians(az1), alt1, c=alt1, cmap='plasma', s=20)
        axs[i, 0].set_title(f"Spiral Grid (n={n})")
        axs[i, 0].set_ylim(90, 0)
        axs[i, 0].set_yticks([0, 30, 60, 90])
        axs[i, 0].set_yticklabels([0, 30, 60, 90])
        cbar1 = plt.colorbar(axs[i, 0].collections[0], ax=axs[i, 0], pad=0.1, shrink=0.8)
        cbar1.set_label("Altitude (deg)")

        # Biased Fibonacci
        r2, th2 = biased_fibonacci_grid(n, bias=bias)
        alt2, az2 = map_to_alt_az(r2, th2, alt_min, alt_max)
        axs[i, 1].scatter(np.radians(az2), alt2, c=alt2, cmap='plasma', s=20)
        axs[i, 1].set_title(f"Biased Fibonacci (n={n})")
        axs[i, 1].set_ylim(90, 0)
        axs[i, 1].set_yticks([0, 30, 60, 90])
        axs[i, 1].set_yticklabels([0, 30, 60, 90])
        cbar2 = plt.colorbar(axs[i, 1].collections[0], ax=axs[i, 1], pad=0.1, shrink=0.8)
        cbar2.set_label("Altitude (deg)")

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    save_path = r"C:\Users\Albert Einstein\Desktop\Einstein docs\Master 2 astro\Grids\spiral_vs_biased_fibonacci.png"
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

# Run
if __name__ == "__main__":
    plot_grids_with_altitude()
