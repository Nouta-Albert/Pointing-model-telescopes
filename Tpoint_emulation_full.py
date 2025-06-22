# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 16:41:38 2025
@author: Albert Einstein
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
import os
import base64
from termcolor import colored

"""This script essentially takes in pointing data, establishes a pointing model based on desired terms and saves the results in the indicated output directory"""

# Fitting parameters
PARAM_BOUNDS = (-1, 1)  # Parameter bounds in radians
CORRELATION_THRESHOLD = 0.8  # Threshold for reporting parameter correlations

# Plotting parameters
PLOT_FIGSIZE = (20, 15)  # Figure size in inches
PLOT_DPI = 300  # DPI for saved figures
PLOT_STYLE = 'default'  # Matplotlib style to use

# ======================================================================
# FUNCTION DEFINITIONS
# ======================================================================

def load_data(file_path):
    """Load pointing data (az_real, el_real, az_actual, el_actual) from a CSV or TXT file."""
    data = pd.read_csv(file_path, sep=r'\s+', header=None, skiprows=0, 
                   names=["az_real", "el_real", "az_actual", "el_actual"])
    # Convert degrees to radians
    data = np.radians(data)
    return data

def angular_difference(a, b):
    """Compute the smallest difference between two angles (in radians)."""
    diff = (a - b + np.pi) % (2 * np.pi) - np.pi
    return diff

def pointing_model(params, az_real, el_real, selected_terms):
    """
    Calculation of the pointing correction, code translated from C.
    Beware the azimuth has a 180 degree offset to standard notations.
    Angles are in radian at input and in arcsec at the output.
    
    In the telescope code it is applied in the following way:
    Azi_final = AZI - Pointing correction azi
    Ele_final = ELE - Pointing correction ele
    
    pointing corrections are the output of this routine.
    """
    # Create a dictionary to map terms to their corresponding parameters
    term_to_param = dict(zip(selected_terms, params))
    
    # Calculate intermediate values
    Zenithal_Angle = 0.5*np.pi - el_real
    COS_Azimut_Angle = np.cos(az_real)
    COS_Elevation_Angle = np.cos(el_real)
    SIN_Azimut_Angle = np.sin(az_real)
    SIN_Elevation_Angle = np.sin(el_real)
    TAN_Elevation_Angle = np.tan(el_real)
    COS_Zenithal_Angle = np.cos(Zenithal_Angle)
    COS_2_Zenithal_Angle = np.cos(2.0 * Zenithal_Angle)
    COS_4_Zenithal_Angle = np.cos(4.0 * Zenithal_Angle)
    SIN_Zenithal_Angle = np.sin(Zenithal_Angle)
    SIN_2_Zenithal_Angle = np.sin(2.0 * Zenithal_Angle)
    SIN_3_Zenithal_Angle = np.sin(3.0 * Zenithal_Angle)
    SIN_4_Zenithal_Angle = np.sin(4.0 * Zenithal_Angle)
    POW_3_Zenithal_Angle = Zenithal_Angle ** 3
    POW_5_Zenithal_Angle = Zenithal_Angle ** 5

    # Initialize corrections
    deltaA = 0.0
    deltaE = 0.0

    # Azimuth corrections
    if "IA" in selected_terms:
        deltaA += term_to_param["IA"]
    if "NPAE" in selected_terms:
        deltaA += term_to_param["NPAE"] * TAN_Elevation_Angle
    if "CA" in selected_terms:
        deltaA += term_to_param["CA"] / COS_Elevation_Angle
    if "AN" in selected_terms:
        deltaA += term_to_param["AN"] * SIN_Azimut_Angle * TAN_Elevation_Angle
    if "AW" in selected_terms:
        deltaA += term_to_param["AW"] * COS_Azimut_Angle * TAN_Elevation_Angle
    if "HSCA" in selected_terms:
        deltaA -= term_to_param["HSCA"] * COS_Azimut_Angle / COS_Elevation_Angle
    if "HVCA" in selected_terms:
        deltaA -= term_to_param["HVCA"] * COS_Azimut_Angle * TAN_Elevation_Angle
    if "ACEC" in selected_terms:
        deltaA -= term_to_param["ACEC"] * COS_Azimut_Angle
    if "ACES" in selected_terms:
        deltaA += term_to_param["ACES"] * SIN_Azimut_Angle
    if "HACA1" in selected_terms:
        deltaA += term_to_param["HACA1"] * COS_Azimut_Angle
    if "HSCZ1" in selected_terms:
        deltaA += term_to_param["HSCZ1"] * COS_Zenithal_Angle / COS_Elevation_Angle
    if "HSCZ2" in selected_terms:
        deltaA += term_to_param["HSCZ2"] * COS_2_Zenithal_Angle / COS_Elevation_Angle
    if "HSSA1" in selected_terms:
        deltaA += term_to_param["HSSA1"] * SIN_Azimut_Angle / COS_Elevation_Angle
    if "HSSZ1" in selected_terms:
        deltaA += term_to_param["HSSZ1"] * SIN_Zenithal_Angle / COS_Elevation_Angle
    if "HSSZ2" in selected_terms:
        deltaA += term_to_param["HSSZ2"] * SIN_2_Zenithal_Angle / COS_Elevation_Angle
    if "HSSZ3" in selected_terms:
        deltaA += term_to_param["HSSZ3"] * SIN_3_Zenithal_Angle / COS_Elevation_Angle
    if "HSCZ4" in selected_terms:
        deltaA += term_to_param["HSCZ4"] * COS_4_Zenithal_Angle / COS_Elevation_Angle
    if "NRX" in selected_terms:
        deltaA += term_to_param["NRX"]
    if "NRY" in selected_terms:
        deltaA += term_to_param["NRY"] * TAN_Elevation_Angle

    # Elevation corrections
    if "IE" in selected_terms:
        deltaE -= term_to_param["IE"]
    if "TF" in selected_terms:
        deltaE += term_to_param["TF"] * COS_Elevation_Angle
    if "TX" in selected_terms:
        deltaE += term_to_param["TX"] / TAN_Elevation_Angle
    if "AN" in selected_terms:
        deltaE += term_to_param["AN"] * COS_Azimut_Angle
    if "AW" in selected_terms:
        deltaE -= term_to_param["AW"] * SIN_Azimut_Angle
    if "ECEC" in selected_terms:
        deltaE += term_to_param["ECEC"] * COS_Elevation_Angle#
    if "ECES" in selected_terms:
        deltaE -= term_to_param["ECES"] * SIN_Elevation_Angle
    if "HZCZ1" in selected_terms:
        deltaE += term_to_param["HZCZ1"] * COS_Zenithal_Angle#
    if "HZCZ2" in selected_terms:
        deltaE += term_to_param["HZCZ2"] * COS_2_Zenithal_Angle
    if "HZCZ4" in selected_terms:
        deltaE += term_to_param["HZCZ4"] * COS_4_Zenithal_Angle
    if "HZSA" in selected_terms:
        deltaE += term_to_param["HZSA"] * SIN_Azimut_Angle#
    if "HZSZ1" in selected_terms:
        deltaE += term_to_param["HZSZ1"] * SIN_Zenithal_Angle#
    if "HZSZ2" in selected_terms:
        deltaE += term_to_param["HZSZ2"] * SIN_2_Zenithal_Angle#
    if "HZSZ4" in selected_terms:
        deltaE += term_to_param["HZSZ4"] * SIN_4_Zenithal_Angle#
    if "PZZ3" in selected_terms:
        deltaE += term_to_param["PZZ3"] * POW_3_Zenithal_Angle
    if "PZZ5" in selected_terms:
        deltaE += term_to_param["PZZ5"] * POW_5_Zenithal_Angle
    if "NRX" in selected_terms:
        deltaE -= term_to_param["NRX"] * SIN_Elevation_Angle
    if "NRY" in selected_terms:
        deltaE += term_to_param["NRY"] * COS_Elevation_Angle

    return deltaA, deltaE

def residuals(params, az_real, el_real, az_actual, el_actual, selected_terms):
    """Compute residuals between actual and modeled positions."""
    deltaA, deltaE = pointing_model(params, az_real, el_real, selected_terms)
    res_az = angular_difference(az_real + deltaA, az_actual)*np.cos(el_real)
    res_el = (el_real + deltaE) - el_actual
    
    return np.concatenate((res_az, res_el))

def calculate_rms(data, fitted_params, selected_terms):
    """Calculate the RMS of each term's contribution and the overall RMS."""
    az_real = data.iloc[:, 0]
    el_real = data.iloc[:, 1]
    az_actual = data.iloc[:, 2]
    el_actual = data.iloc[:, 3]
    
    # Calculate overall residuals
    deltaA, deltaE = pointing_model(fitted_params, az_real, el_real, selected_terms)
    res_az = angular_difference(az_real + deltaA, az_actual)*np.cos(el_real)
    res_el = (el_real + deltaE) - el_actual
    overall_rms_az = np.sqrt(np.mean(res_az**2)) * (180/np.pi) * 3600  # Convert to arcseconds
    overall_rms_el = np.sqrt(np.mean(res_el**2)) * (180/np.pi) * 3600
    
    return overall_rms_az, overall_rms_el

def fit_pointing_model(data, selected_terms, initial_params):
    """Perform least-squares optimization and return parameters + uncertainties."""
    az_real = data.iloc[:, 0]
    el_real = data.iloc[:, 1]
    az_actual = data.iloc[:, 2]
    el_actual = data.iloc[:, 3]
    
    result = least_squares(
        residuals, 
        initial_params, 
        args=(az_real, el_real, az_actual, el_actual, selected_terms), 
        method='trf', 
        bounds=PARAM_BOUNDS
    )
    
    # Compute uncertainties and covariance
    mse = np.sum(result.fun**2) / (0.5*len(result.fun) - len(result.x))
    cov = mse * np.linalg.inv(result.jac.T @ result.jac)
    std_errors = np.degrees(np.sqrt(np.diag(cov)))*3600
    
    return result, std_errors, cov

def analyze_correlations(cov_matrix, selected_terms):
    """
    Compute correlations and return strongly correlated pairs.
    """
    # Compute correlation matrix
    std_devs = np.sqrt(np.diag(cov_matrix))
    corr_matrix = cov_matrix / np.outer(std_devs, std_devs)
    
    # Find strongly correlated pairs
    strong_corrs = []
    for i in range(len(selected_terms)):
        for j in range(i+1, len(selected_terms)):
            if abs(corr_matrix[i,j]) > CORRELATION_THRESHOLD:
                strong_corrs.append(
                    (selected_terms[i], selected_terms[j], corr_matrix[i,j])
                )
    
    return corr_matrix, strong_corrs

def analyze_and_plot_all(data, fitted_params, selected_terms, output_dir, cov_matrix):
    """
    Unified plotting function with 3x3 grid layout.
    """
    # --- Data Preparation ---
    az_real = data.iloc[:, 0]
    el_real = data.iloc[:, 1]
    az_actual = data.iloc[:, 2]
    el_actual = data.iloc[:, 3]
    
    # Compute residuals
    deltaA, deltaE = pointing_model(fitted_params, az_real, el_real, selected_terms)
    res_az = angular_difference(az_real + deltaA, az_actual) * np.cos(el_real)
    res_el = (el_real + deltaE) - el_actual
    
    res_az_arcsec = res_az * (180/np.pi) * 3600
    res_el_arcsec = res_el * (180/np.pi) * 3600
    rms_az = np.sqrt(np.mean(res_az_arcsec**2))
    rms_el = np.sqrt(np.mean(res_el_arcsec**2))
    rms_tot = np.sqrt(rms_az**2 + rms_el**2)
    
    # Compute correlation matrix
    std_devs = np.sqrt(np.diag(cov_matrix))
    corr_matrix = cov_matrix / np.outer(std_devs, std_devs)

    # --- Create Figure ---
    plt.style.use(PLOT_STYLE)
    fig = plt.figure(figsize=PLOT_FIGSIZE)
    fig.suptitle(f'Pointing Model Analysis: RMS ΔAz: {rms_az:.2f}"  RMS ΔEl: {rms_el:.2f}"  RMS  Δtot: {rms_tot:.2f}"', 
                fontsize=15, y=0.98)
    
    # 3x3 Grid
    gs = fig.add_gridspec(3, 3, hspace=0.4, wspace=0.3)

    # --- Row 1: Polar Plots + Heatmap --
    # Polar Plot 1: Point distribution
    ax1 = fig.add_subplot(gs[0, 0], projection='polar')
    r = np.degrees(el_real)   # Altitude
    theta = az_real           # Azimuth (radians)
    
    sc = ax1.scatter(theta, r, s=10, c=r, cmap='cool', alpha=0.7)
    plt.colorbar(sc, ax=ax1, label='Altitude (degrees)', pad=0.1)
    
    ax1.set_theta_zero_location('N')
    ax1.set_theta_direction(-1)  # Counterclockwise
    ax1.set_ylim(90, 0)         # Zenith at center
    ax1.set_yticks([90, 60, 30, 0])
    ax1.set_title("Points Distribution", pad=15)
    ax1.grid(alpha=0.3)

    # Polar Plot 2: Residuals
    ax2 = fig.add_subplot(gs[0, 1])
    residual_mag = np.sqrt(res_az_arcsec**2 + res_el_arcsec**2)
    circle_80 = np.percentile(residual_mag, 80)
    circle_100 = np.percentile(residual_mag, 100)
    
    ax2.scatter(res_az_arcsec, res_el_arcsec, alpha=0.7)
    ax2.axhline(y=0, color='r', linestyle='--', linewidth=1)
    ax2.axvline(x=0, color='r', linestyle='--', linewidth=1)
    
    circle_80_artist = plt.Circle((0, 0), circle_80, fill=False, color='green', linestyle='--', alpha=0.7)
    circle_100_artist = plt.Circle((0, 0), circle_100, fill=False, color='red', linestyle='--', alpha=0.7)
    ax2.add_artist(circle_80_artist)
    ax2.add_artist(circle_100_artist)
    
    ax2.annotate(f'{circle_80:.1f}"', xy=(circle_80, 0), xytext=(5, 5), 
                 textcoords='offset points', color='green', fontsize=10)
    ax2.annotate(f'{circle_100:.1f}"', xy=(circle_100, 0), xytext=(5, 5), 
                 textcoords='offset points', color='red', fontsize=10)
    ax2.set_xlabel('ΔAz (arcsec)')
    ax2.set_ylabel('ΔEl (arcsec)')
    ax2.set_title('Residuals: ΔAz vs ΔEl')
    ax2.grid(True, linestyle=':', alpha=0.7)
    
    max_limit = max(abs(res_az_arcsec).max(), abs(res_el_arcsec).max()) * 1.1
    ax2.set_xlim(-max_limit, max_limit)
    ax2.set_ylim(-max_limit, max_limit)
    ax2.set_aspect('equal')

    # Heatmap
    ax3 = fig.add_subplot(gs[0, 2])
    im = ax3.imshow(corr_matrix, cmap='coolwarm', vmin=-1, vmax=1)
    plt.colorbar(im, ax=ax3, label='Correlation')
    ax3.set_xticks(range(len(selected_terms)))
    ax3.set_yticks(range(len(selected_terms)))
    ax3.set_xticklabels(selected_terms, rotation=90)
    ax3.set_yticklabels(selected_terms)
    ax3.set_title("Parameter Correlations", pad=15)

    # --- Row 2: ΔAz Diagnostics ---
    ax4 = fig.add_subplot(gs[1, 0])
    ax4.scatter(np.degrees(az_real), res_az_arcsec, s=10, color='blue', alpha=0.7)
    ax4.set_xlabel("Azimuth [deg]"); ax4.set_ylabel("ΔAz [arcsec]")
    ax4.grid(alpha=0.3)

    ax5 = fig.add_subplot(gs[1, 1])
    ax5.scatter(np.degrees(el_real), res_az_arcsec, s=10, color='green', alpha=0.7)
    ax5.set_xlabel("Elevation [deg]"); ax5.set_ylabel("ΔAz [arcsec]")
    ax5.grid(alpha=0.3)

    ax6 = fig.add_subplot(gs[1, 2])
    ax6.hist(res_az_arcsec, bins=30, color='blue', alpha=0.7, density=True)
    ax6.set_xlabel("ΔAz [arcsec]"); ax6.grid(alpha=0.3)

    # --- Row 3: ΔEl Diagnostics ---
    ax7 = fig.add_subplot(gs[2, 0])
    ax7.scatter(np.degrees(az_real), res_el_arcsec, s=10, color='red', alpha=0.7)
    ax7.set_xlabel("Azimuth [deg]"); ax7.set_ylabel("ΔEl [arcsec]")
    ax7.grid(alpha=0.3)

    ax8 = fig.add_subplot(gs[2, 1])
    ax8.scatter(np.degrees(el_real), res_el_arcsec, s=10, color='purple', alpha=0.7)
    ax8.set_xlabel("Elevation [deg]"); ax8.set_ylabel("ΔEl [arcsec]")
    ax8.grid(alpha=0.3)

    ax9 = fig.add_subplot(gs[2, 2])
    ax9.hist(res_el_arcsec, bins=30, color='red', alpha=0.7, density=True)
    ax9.set_xlabel("ΔEl [arcsec]"); ax9.grid(alpha=0.3)

    # Save figure
    plt.savefig(os.path.join(output_dir, "Results plots.png"), 
               dpi=PLOT_DPI, bbox_inches='tight')
    plt.close()

def save_results(output_dir, selected_terms, result, std_errors, overall_rms_az, overall_rms_el, corr_matrix, strong_corrs):
    """Save all results to files."""
    # Save parameters and RMS values
    with open(os.path.join(output_dir, "pointing_results.txt"), "w", encoding='utf-8') as f:
        f.write("Fitted Pointing Terms with 1-sigma uncertainties (arcsec):\n")
        for term, value, error in zip(selected_terms, result.x, std_errors):
            f.write(f"{term}: {np.degrees(value) * 3600:.4f} ± {error:.4f}\n")
        f.write(f"\nFit RMS Az: {overall_rms_az:.2f} arcsec\n")
        f.write(f"Fit RMS El: {overall_rms_el:.2f} arcsec\n")
        f.write(f"Fit total RMS: {np.sqrt(overall_rms_az**2 + overall_rms_el**2):.2f} arcsec\n")
        pop_sd = np.sqrt(np.sum(result.fun**2)/(0.5*len(result.fun) - len(result.x))) * (180/np.pi) * 3600
        f.write(f"\nPopulation standard deviation: {pop_sd:.2f} arcsec\n")
        # Save correlation information if any strong correlations exist
        if strong_corrs:
            f.write("\nStrongly correlated parameters (|r| > 0.8):\n")
            for term1, term2, corr in strong_corrs:
                f.write(f"{term1} vs {term2}: {corr:.2f}\n")
    
    # Save full correlation matrix
    with open(os.path.join(output_dir, "correlation_matrix.csv"), "w") as f:
        header = "Parameter".ljust(12) + "," + ",".join([f"{term:^10}" for term in selected_terms])
        f.write(header + "\n")
        f.write("-" * (12 + 11 * len(selected_terms)) + "\n")
        for i, term in enumerate(selected_terms):
            row_str = f"{term.ljust(12)},"
            corr_values = [f"{corr_matrix[i,j]:^10.3f}" for j in range(len(selected_terms))]
            row_str += ",".join(corr_values)
            f.write(row_str + "\n")
        f.write("\nNotes:\n")
        f.write("Diagonal values (self-correlations) are always 1.000\n")
        f.write(f"Values > |{CORRELATION_THRESHOLD:.1f}| may indicate significant parameter correlations\n")
    #combine results into one html file
    bundle_to_html(output_dir)

def bundle_to_html(output_dir):
    """Combine existing files into one self-contained HTML file"""

    # Read and encode image as Base64
    img_path = os.path.join(output_dir, "Results plots.png")
    with open(img_path, "rb") as img_file:
        encoded_img = base64.b64encode(img_file.read()).decode("utf-8")
    img_data_uri = f"data:image/png;base64,{encoded_img}"

    # Read text files
    with open(os.path.join(output_dir, "pointing_results.txt"), 'r') as f:
        pointing_results = f.read()

    with open(os.path.join(output_dir, "correlation_matrix.csv"), 'r') as f:
        correlation_matrix = f.read()

    # HTML content
    html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Pointing Model Results</title>
    <style>
        body {{ font-family: Arial; margin: 20px }}
        .section {{ margin-bottom: 30px; border-bottom: 1px solid #eee; padding-bottom: 15px }}
        img {{ max-width: 100% }}
    </style>
</head>
<body>
    <h1>Pointing Model Results</h1>

    <div class="section">
        <h2>Diagnostic Plots</h2>
        <img src="{img_data_uri}">
    </div>

    <div class="section">
        <h2>Pointing Results</h2>
        <pre>{pointing_results}</pre>
    </div>

    <div class="section">
        <h2>Correlation Matrix</h2>
        <pre>{correlation_matrix}</pre>
    </div>
</body>
</html>
    """

    # Save combined HTML
    with open(os.path.join(output_dir, "combined_results.html"), 'w') as f:
        f.write(html_content)

    # Delete temporary files
    for file in ["Results plots.png", "pointing_results.txt", "correlation_matrix.csv"]:
        os.remove(os.path.join(output_dir, file))


def main(input_file, output_dir, selected_terms):
    os.makedirs(output_dir, exist_ok=True)
    data = load_data(input_file)
    print("Fitting the following terms:", selected_terms)

    initial_params = np.zeros(len(selected_terms))
    result, std_errors, cov = fit_pointing_model(data, selected_terms, initial_params)

    corr_matrix, strong_corrs = analyze_correlations(cov, selected_terms)
    overall_rms_az, overall_rms_el = calculate_rms(data, result.x, selected_terms)

    print(f"\nFit RMS Az: {overall_rms_az:.2f} arcsec")
    print(f"Fit RMS El: {overall_rms_el:.2f} arcsec")
    print(f"Fit total RMS: {np.sqrt(overall_rms_az**2 + overall_rms_el**2):.2f} arcsec")

    pop_std = np.sqrt(np.sum(result.fun**2)/(0.5*len(result.fun) - len(result.x))) * (180/np.pi) * 3600
    print(f"\nPopulation Standard Deviation: {pop_std:.2f} arcsec\n")

    for term, value, error in zip(selected_terms, result.x, std_errors):
        print(colored(f"{term}: {np.degrees(value) * 3600:.4f} ± {error:.4f}", "magenta"))

    if strong_corrs:
        print("\nStrongly correlated parameters (|r| > 0.8):")
        for term1, term2, corr in strong_corrs:
            print(f"{term1} vs {term2}: {corr:.2f}")
    else:
        print("\nNo strongly correlated parameter pairs found (all |r| ≤ 0.8)")

    analyze_and_plot_all(data, result.x, selected_terms, output_dir, cov)
    save_results(output_dir, selected_terms, result, std_errors, overall_rms_az, overall_rms_el, corr_matrix, strong_corrs)
    

    print("\nResults saved in:", output_dir)


def run_pointing_model(input_file, output_dir, selected_terms):
    try:
        main(input_file, output_dir, selected_terms)
        print(colored("\nPointing model successfully computed and saved. Check output directory for results!", "green"))
    except Exception as e:
        print(colored(f"Failed to compute pointing model: {e}", "red"))
