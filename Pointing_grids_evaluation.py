# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 12:23:46 2025

@author: Albert Einstein
"""

# -*- coding: utf-8 -*-
"""
Enhanced Pointing Model Grid Evaluation - Fixed Version
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
from tqdm import tqdm

# Import your functions
from Tpoint_emulation_full import (
    pointing_model, fit_pointing_model,
    calculate_rms, analyze_and_plot_all
)
from Grids_comparison import (
    spiral_grid, biased_fibonacci_grid,
    map_to_alt_az
)

# --- Constants ---
TRUE_PARAMS = np.array([5.3268, 42.9147, 95.6802, 291.2327, 61.5283,
                        13.0052, 81.6729, -133.4021, 29.4735]) * np.pi / (180 * 3600)
SELECTED_TERMS = ["AN", "AW", "IA", "IE", "CA", "TF", "PZZ3", "PZZ5", "HZCZ4"]
OUTPUT_DIR = r"C:\Users\Albert Einstein\Desktop\Einstein docs\Master 2 astro\Grids"

def generate_simulated_data(n_points, grid_type, noise_level=2, add_noise=True):
    """Generate simulated observations with proper noise modeling"""
    if grid_type == 'spiral':
        r, th = spiral_grid(n_points)
    else:
        r, th = biased_fibonacci_grid(n_points, bias=0.8)
    
    el_real, az_real = map_to_alt_az(r, th)  # in degrees
    az_rad, el_rad = np.radians(az_real), np.radians(el_real)
    
    # Apply true pointing model
    deltaA, deltaE = pointing_model(TRUE_PARAMS, az_rad, el_rad, SELECTED_TERMS)
    
    # Add Gaussian noise (more realistic than uniform)
    if add_noise:
        noise_az = np.random.normal(0, noise_level/3600, n_points) / np.cos(el_rad)
        noise_el = np.random.normal(0, noise_level/3600, n_points)
    else:
        noise_az = noise_el = np.zeros(n_points)
    
    return pd.DataFrame({
        'az_real': az_real,
        'el_real': el_real,
        'az_actual': np.degrees(az_rad + deltaA) + noise_az,
        'el_actual': np.degrees(el_rad + deltaE) + noise_el
    })

def run_simulation(n_points, grid_type, n_trials=2000):
    """Run multiple trials for statistical significance"""
    results = {
        'rms_az': [],
        'rms_el': [],
        'params': [],
        'std_errors': []
    }
    
    for _ in tqdm(range(n_trials), desc=f"{grid_type} {n_points}pts"):
        data_deg = generate_simulated_data(n_points, grid_type)
        data_rad = data_deg.apply(np.radians)
        result, std_errors, _ = fit_pointing_model(data_rad, SELECTED_TERMS)
        rms_az, rms_el = calculate_rms(data_rad, result.x, SELECTED_TERMS)
        
        results['rms_az'].append(rms_az)
        results['rms_el'].append(rms_el)
        results['params'].append(result.x)
        results['std_errors'].append(std_errors)
    
    # Calculate total RMS for each trial
    results['total_rms'] = np.sqrt(np.array(results['rms_az'])**2 + np.array(results['rms_el'])**2)
    
    return results

def plot_rms_comparison(all_results):
    """Visualize RMS results across configurations"""
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    
    point_counts = [20, 60, 100]
    metrics = ['rms_az', 'rms_el', 'total_rms']
    titles = ['Azimuth RMS (arcsec)', 'Elevation RMS (arcsec)', 'Total RMS (arcsec)']
    
    for i, metric in enumerate(metrics):
        for j, n in enumerate(point_counts):
            spiral_data = all_results['spiral'][n][metric]
            fib_data = all_results['fibonacci'][n][metric]
            
            positions = [j*3, j*3+1]
            bp1 = ax[i].boxplot([spiral_data], positions=[positions[0]], widths=0.6,
                               patch_artist=True,
                               boxprops=dict(facecolor='blue', alpha=0.5),
                               whiskerprops=dict(color='blue'),
                               capprops=dict(color='blue'),
                               medianprops=dict(color='darkblue'),
                               flierprops=dict(marker='o', markersize=5, alpha=0.5))
            
            bp2 = ax[i].boxplot([fib_data], positions=[positions[1]], widths=0.6,
                               patch_artist=True,
                               boxprops=dict(facecolor='orange', alpha=0.5),
                               whiskerprops=dict(color='orange'),
                               capprops=dict(color='orange'),
                               medianprops=dict(color='darkorange'),
                               flierprops=dict(marker='o', markersize=5, alpha=0.5))
        
        ax[i].set_title(titles[i])
        ax[i].set_xticks([0.5, 3.5, 6.5], point_counts)
        ax[i].set_xlabel('Number of Points')
        ax[i].set_ylabel('RMS (arcsec)')
        ax[i].legend([bp1["boxes"][0], bp2["boxes"][0]], ['Spiral', 'Fibonacci'])
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'rms_comparison.png'), dpi=150)
    plt.close()
    
def plot_pointing_corrections():
    """Create 2D maps of pointing corrections"""
    az = np.linspace(0, 2*np.pi, 100)  # 0-360° in radians
    el = np.linspace(np.radians(30), np.radians(85), 100)  # 30-85° elevation
    az_grid, el_grid = np.meshgrid(az, el)
    
    # Calculate corrections (vectorized)
    deltaA, deltaE = pointing_model(TRUE_PARAMS, az_grid.flatten(), el_grid.flatten(), SELECTED_TERMS)
    deltaA = deltaA.reshape(az_grid.shape) * (180/np.pi) * 3600  # arcsec
    deltaE = deltaE.reshape(el_grid.shape) * (180/np.pi) * 3600
    delta_total = np.sqrt(deltaA**2 + deltaE**2)
    
    # Create polar plots
    fig = plt.figure(figsize=(18, 6))
    titles = ['ΔAz Correction (arcsec)', 'ΔEl Correction (arcsec)', 'Total Correction (arcsec)']
    data = [deltaA, deltaE, delta_total]
    
    for i in range(3):
        ax = fig.add_subplot(1, 3, i+1, projection='polar')
        ax.set_theta_zero_location('N')  # 0° at North
        ax.set_theta_direction(-1)       # Clockwise azimuth
        
        # Convert elevation to zenith distance (90-el) for proper polar display
        r = 90 - np.degrees(el_grid)
        
        # Plot with color mapping
        pc = ax.pcolormesh(az_grid, r, data[i], 
                          cmap='plasma', 
                          shading='auto',
                          vmin=np.min(data[i]), 
                          vmax=np.max(data[i]))
        
        plt.colorbar(pc, ax=ax, label='arcsec', pad=0.1)
        ax.set_title(titles[i], pad=20)
        ax.set_ylim(0, 60)  # 30-90° elevation
        ax.set_yticks([0, 30, 60])
        ax.set_yticklabels(['90°', '60°', '30°'])  # Zenith at center
        
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'polar_corrections.png'), dpi=150)
    plt.close()

def save_statistical_results(all_results):
    """Save comprehensive results to CSV"""
    for grid_type in ['spiral', 'fibonacci']:
        for n_points in [20, 60, 100]:
            results = all_results[grid_type][n_points]
            
            # Calculate statistics
            mean_params = np.mean(results['params'], axis=0)
            std_params = np.std(results['params'], axis=0)
            mean_std_errors = np.mean(results['std_errors'], axis=0)
            
            # Create DataFrame
            df = pd.DataFrame({
                'Term': SELECTED_TERMS,
                'True_Value': TRUE_PARAMS * (180/np.pi) * 3600,
                'Mean_Fitted': mean_params * (180/np.pi) * 3600,
                'StdDev_Fitted': std_params * (180/np.pi) * 3600,
                'Mean_Uncertainty': mean_std_errors
            })
            
            # Add RMS stats
            rms_stats = pd.DataFrame({
                'Metric': ['Az RMS', 'El RMS', 'Total RMS'],
                'Mean': [
                    np.mean(results['rms_az']),
                    np.mean(results['rms_el']),
                    np.mean(results['total_rms'])
                ],
                'StdDev': [
                    np.std(results['rms_az']),
                    np.std(results['rms_el']),
                    np.std(results['total_rms'])
                ]
            })
            
            # Save to CSV
            dir_path = os.path.join(OUTPUT_DIR, f"{grid_type}_{n_points}pts")
            os.makedirs(dir_path, exist_ok=True)
            
            df.to_csv(os.path.join(dir_path, 'parameters.csv'), index=False)
            rms_stats.to_csv(os.path.join(dir_path, 'rms_stats.csv'), index=False)

def main():
    """Run comprehensive analysis"""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Run simulations for all configurations
    all_results = {'spiral': {}, 'fibonacci': {}}
    
    for n_points in [20, 60, 100]:
        print(f"\n=== Running simulations for {n_points} points ===")
        all_results['spiral'][n_points] = run_simulation(n_points, 'spiral')
        all_results['fibonacci'][n_points] = run_simulation(n_points, 'fibonacci')
    
    # Generate plots and save results
    plot_rms_comparison(all_results)
    plot_pointing_corrections()
    save_statistical_results(all_results)
    
    print("\n=== Analysis Complete ===")
    print(f"Results saved to: {OUTPUT_DIR}")
    print("Includes:")
    print("- RMS comparison plots")
    print("- Pointing correction maps")
    print("- Statistical results for each configuration")

if __name__ == "__main__":
    main()