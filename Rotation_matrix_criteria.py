import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# === USER INPUT ===
directory = "/home/astro/alberte2/PM_test_LPM/"
output_plot_name = "/home/astro/alberte2/PM_test_LPM/cd_matrix_validation.png"
output_text_name = "/home/astro/alberte2/PM_test_LPM/cd_matrix_validation.txt"

# === FUNCTIONS ===

def polar_decompose_cd(cd):
    """Decompose CD into rotation (R) and symmetric stretch/shear (P)."""
    U, s, Vt = np.linalg.svd(cd)
    R = U @ Vt
    P = Vt.T @ np.diag(s) @ Vt
    return R, P

def compute_angle_from_matrix(matrix):
    """Compute rotation angle (degrees) from 2x2 matrix."""
    theta_rad = np.arctan2(-matrix[0, 1], matrix[1, 1])
    return np.degrees(theta_rad)

def analyze_cd_matrix(cd):
    """Perform full analysis: naive and true angles, delta, and distortion."""
    angle_naive = compute_angle_from_matrix(cd)
    R, P = polar_decompose_cd(cd)
    angle_true = compute_angle_from_matrix(R)

    delta_angle = np.abs(angle_true - angle_naive)

    # Measure anisotropy from P
    eigvals = np.linalg.eigvals(P)
    eccentricity = np.max(eigvals) / np.min(eigvals)

    return angle_naive, angle_true, delta_angle, eccentricity

def process_directory(directory):
    filenames = []
    naive_angles = []
    true_angles = []
    deltas = []
    eccentricities = []

    for fname in sorted(os.listdir(directory)):
        if not fname.lower().endswith(".solved.fits"):
            continue

        fpath = os.path.join(directory, fname)
        try:
            with fits.open(fpath) as hdul:
                header = hdul[0].header
                cd = np.array([
                    [header['CD1_1'], header['CD1_2']],
                    [header['CD2_1'], header['CD2_2']]
                ])

                angle_naive, angle_true, delta_angle, ecc = analyze_cd_matrix(cd)

                filenames.append(fname)
                naive_angles.append(angle_naive)
                true_angles.append(angle_true)
                deltas.append(delta_angle)
                eccentricities.append(ecc)

        except Exception as e:
            print(f"Skipping {fname}: {e}")

    return filenames, naive_angles, true_angles, deltas, eccentricities

def save_summary(filenames, naive_angles, true_angles, deltas, eccentricities, path):
    with open(path, "w") as f:
        f.write(f"{'Filename':<40} {'θ_naive(deg)':>12} {'θ_true(deg)':>12} {'Δθ(deg)':>10} {'Eccentricity':>15}\n")
        f.write("-" * 95 + "\n")
        for fname, a1, a2, dtheta, ecc in zip(filenames, naive_angles, true_angles, deltas, eccentricities):
            f.write(f"{fname:<40} {a1:12.3f} {a2:12.3f} {dtheta:10.3f} {ecc:15.6f}\n")
    print(f"[✔] Summary saved to: {path}")

def plot_metrics(filenames, deltas, eccentricities, output_path):
    x = np.arange(len(filenames))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    ax1.plot(x, deltas, marker='o', linestyle='-', color='tab:blue')
    ax1.set_title("Rotation Angle Error (Δθ)")
    ax1.set_xlabel("File Index")
    ax1.set_ylabel("Angle Difference (deg)")
    ax1.grid(True, linestyle=':')

    ax2.plot(x, eccentricities, marker='s', linestyle='-', color='tab:red')
    ax2.set_title("Distortion Measure (Eccentricity)")
    ax2.set_xlabel("File Index")
    ax2.set_ylabel("λ₁ / λ₂ (Scaling Ratio)")
    ax2.grid(True, linestyle=':')

    fig.suptitle("CD Matrix Validation: Rotation Angle Reliability", fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(output_path, dpi=300)
    print(f"[✔] Plot saved to: {output_path}")
    plt.show()

# === MAIN EXECUTION ===

if __name__ == "__main__":
    filenames, naive, true, delta, ecc = process_directory(directory)

    text_out = os.path.join(directory, output_text_name)
    plot_out = os.path.join(directory, output_plot_name)

    save_summary(filenames, naive, true, delta, ecc, text_out)
    plot_metrics(filenames, delta, ecc, plot_out)
