import os
import subprocess
import pandas as pd
from datetime import datetime
from astropy.io import fits
from astropy.coordinates import EarthLocation, SkyCoord, AltAz
from astropy.time import Time
from astropy.wcs import WCS
import numpy as np
import shutil
import glob
import time
from termcolor import colored
import re

################################################################################
# Directories
#archives_dir = "/data/ECAMRAW/raw/"
archives_dir = "/home/astro/alberte2/euler02_data_ECAM/"
temp_working_dir = "/home/astro/alberte2/Archives_ECAM/"
solved_output_dir = os.path.join(temp_working_dir, "solved/")
os.makedirs(solved_output_dir, exist_ok=True)
output_dir = "/home/astro/alberte2/Results_ECAM_analysis_Feb-May/"
#################################################################################

# Output files
unsolved_files_log = os.path.join(output_dir, "unsolved_files.txt")
solving_stats_file = os.path.join(output_dir, "solving_stats.csv")
pointing_errors_file = os.path.join(output_dir, "pointing_errors.csv")
pointing_data_file = os.path.join(output_dir, "pointing_data.csv")

# Observatory location for calculations
obs_location = EarthLocation(lon=-70.732951945082, lat=-29.259497822801, height=2378.0)

# Initialize log and CSV files
with open(unsolved_files_log, 'w') as f:
    f.write("List of unsolved files:\n")

if not os.path.exists(solving_stats_file):
    pd.DataFrame(columns=['File', 'Solve_Time']).to_csv(solving_stats_file, index=False)

if not os.path.exists(pointing_errors_file):
    pd.DataFrame(columns=['Timestamp', 'az_error', 'el_error']).to_csv(pointing_errors_file, index=False)

if not os.path.exists(pointing_data_file):
    pd.DataFrame(columns=['Obs_date', 'az_expected', 'el_expected', 'az_actual', 'el_actual']).to_csv(pointing_data_file, index=False)

# Function to get directories in chronological order
def get_sorted_directories(base_dir, start_date=None, end_date=None):
    dirs = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    
    if start_date and end_date:
        start_date = datetime.strptime(start_date, "%Y-%m-%d")
        end_date = datetime.strptime(end_date, "%Y-%m-%d")
    
    filtered_dirs = []
    for d in dirs:
        match = re.match(r"(\d{4}-\d{2}-\d{2})", d)
        if match:
            dir_date = datetime.strptime(match.group(1), "%Y-%m-%d")
            if start_date <= dir_date <= end_date:
                filtered_dirs.append(d)
        else:
            print(f"Skipping invalid directory: {d}")
    
    return [os.path.join(base_dir, d) for d in sorted(filtered_dirs, reverse=True)]

# Function to filter and identify science images
def is_science_image(fits_file):
    try:
        with fits.open(fits_file) as hdul:
            header = hdul[0].header
            image = header.get('OBJECT', '')
            exposure = header.get('HIERARCH OGE OBS TYPE', '')
            readout_port = header.get('HIERARCH OGE DET OUT RNAME', '')
            return (image not in ['', 'bias', 'flat', 'dark', 'led'] and exposure != 'ABTR' and readout_port != 'ALL')
    except Exception as e:
        print(f"Error reading FITS file {fits_file}: {e}")
        return False

# Function to perform platesolving and collect results
def plate_solve_and_log(fits_file):
    file_name = os.path.basename(fits_file)
    output_name = file_name.replace(".fits", ".solved.fits")
    solved_file_path = os.path.join(solved_output_dir, output_name)
    
    with fits.open(fits_file) as hdul:
        header = hdul[0].header
        ra_approx = header.get('RA', None)
        dec_approx = header.get('DEC', None)
        
        if ra_approx is None or dec_approx is None:
            print(colored(f"Skipping {file_name}: Missing RA/Dec in header.", "yellow"))
            return False
        
    try:
        start_time = time.time()
        result = subprocess.run([
            "solve-field "
            "--overwrite "
            "--no-verify "
            f"--ra {ra_approx} "
            f"--dec {dec_approx} "
            "--scale-units arcsecperpix "
            "--scale-low 0.2 "
            "--scale-high 0.3 "
            "--radius 3 "
            "--depth 100 "
            "--objs 200 "
            "--no-plots "
            f"--dir {solved_output_dir} "
            f"--new-fits {output_name} {fits_file}"
        ], cwd=solved_output_dir, shell=True, capture_output=True, text=True, check=True)
        
        end_time = time.time()
        elapsed_time = end_time - start_time
        #print(colored(f"Solved {file_name} in {elapsed_time:.2f} seconds.", "green"))

        # Append solving stats
        # Check if solved file was created
        if os.path.exists(solved_file_path):
            print(colored(f"Solved {file_name} in {elapsed_time:.2f} seconds.", "green"))
            pd.DataFrame([{ 'File': fits_file, 'Solve_Time': elapsed_time }]).to_csv(solving_stats_file, mode='a', header=False, index=False)
            return True
        else:
            print(colored(f"No .solved.fits file found for {file_name}. Skipping logging.", "red"))
            with open(unsolved_files_log, 'a') as f:
                f.write(f"{fits_file}\n")
            return False
    except subprocess.CalledProcessError as e:
        print(colored(f"Failed to solve: {file_name}.", "red"))
        with open(unsolved_files_log, 'a') as f:
            f.write(f"{fits_file}\n")
        return False

# Function to calculate pointing errors and save data
def calculate_and_log_pointing_errors():
    pointing_data = []
    pointing_errors = []

    for file_name in os.listdir(solved_output_dir):
        if file_name.endswith('.solved.fits'):
            fits_path = os.path.join(solved_output_dir, file_name)

            try:
                with fits.open(fits_path) as hdul:
                    header = hdul[0].header
                    w = WCS(header)
                    obs_time = Time(float(header.get('MJD-OBS')), format='mjd')
                    obs_date = header.get('DATE-OBS')
                    
                    x_ref = 2121 #header.get('_RPIX1', None)
                    y_ref = 2121.5 #header.get('_RPIX2', None)
                    if x_ref is None or y_ref is None:
                        continue

                    sky_coord = SkyCoord(ra=header.get('RA', None), dec=header.get('DEC', None), frame="icrs", unit="deg")
                    altaz_frame = AltAz(location=obs_location, obstime=obs_time)
                    altaz_coord = sky_coord.transform_to(altaz_frame)
                    az_expected = altaz_coord.az.deg
                    el_expected = altaz_coord.alt.deg

                    pixel_world_coords = w.pixel_to_world(x_ref, y_ref)
                    sky_coord = SkyCoord(ra=pixel_world_coords.ra, dec=pixel_world_coords.dec, frame="icrs", unit="deg")
                    altaz_coord = sky_coord.transform_to(altaz_frame)
                    az_actual = altaz_coord.az.deg
                    el_actual = altaz_coord.alt.deg

                    el_err = (el_actual - el_expected) * 3600
                    az_err = (az_actual - az_expected) * np.cos(np.radians(el_expected)) * 3600

                    pointing_data.append({'Obs_date': obs_date, 'az_expected': az_expected, 'el_expected': el_expected, 'az_actual': az_actual, 'el_actual': el_actual})
                    pointing_errors.append({'Timestamp': obs_date, 'az_error': az_err, 'el_error': el_err})
            except Exception as e:
                print(f"Error processing {file_name}: {e}")

    # Save data
    pd.DataFrame(pointing_data).to_csv(pointing_data_file, mode='a', header=False, index=False)
    pd.DataFrame(pointing_errors).to_csv(pointing_errors_file, mode='a', header=False, index=False)

##############################################################
#                                                            #
#                         Execution                          #
#                                                            #
##############################################################

# Initialize counters
solved_files = 0
total_files = 0
solving_times = []

#########################################################################
# User-defined start and end dates
#start_date = input("Enter start date (YYYY-MM-DD): ")
#end_date = input("Enter end date (YYYY-MM-DD): ")
#start_date = "2010-09-16"
start_date = "2025-01-01"
end_date = "2025-05-12"
#########################################################################

# Loop through directories chronologically
directories = get_sorted_directories(archives_dir, start_date, end_date)

for day_dir in directories:
    print(colored(f"Processing directory: {day_dir}", "cyan"))
    
    # Collect and sort .fits files based on observation date
    fits_files = []
    for root, _, files in os.walk(day_dir):
        for file in files:
            if file.endswith(".fits"):
                file_path = os.path.join(root, file)
                
                # Read observation date from the FITS header
                try:
                    with fits.open(file_path) as hdul:
                        obs_date = hdul[0].header.get('DATE-OBS', '')
                        if obs_date:
                            obs_date_dt = datetime.strptime(obs_date, "%Y-%m-%dT%H:%M:%S.%f")
                            fits_files.append((obs_date_dt, file_path))
                except Exception as e:
                    print(colored(f"Error reading {file}: {e}", "red"))

    # Sort the files by observation date
    fits_files.sort(key=lambda x: x[0], reverse=True)  # Most recent first

    # Process each sorted FITS file
    for _, file_path in fits_files:
        if is_science_image(file_path):
            total_files += 1
            print(colored(f"Platesolving: {file_path}", "yellow"))
            temp_file_path = os.path.join(temp_working_dir, os.path.basename(file_path))
            shutil.copy(file_path, temp_file_path)
            
            if plate_solve_and_log(temp_file_path):
                solved_files += 1
                calculate_and_log_pointing_errors()
            
            solving_times.append(time.time() - os.path.getmtime(temp_file_path))  # Track solving time
            for file in os.listdir(solved_output_dir):
                os.remove(os.path.join(solved_output_dir, file))
            # Cleanup: Delete all files and directories in temp_working_dir except the "solved" folder
            for item in os.listdir(temp_working_dir):
                item_path = os.path.join(temp_working_dir, item)
                if item != "solved":  # Exclude the "solved" folder
                    if os.path.isfile(item_path):
                        os.remove(item_path)
                    elif os.path.isdir(item_path):
                        shutil.rmtree(item_path)



# Calculate and display success statistics
if total_files > 0:
    success_percentage = (solved_files / total_files) * 100
    max_time = np.max(solving_times)
    min_time = np.min(solving_times)
    avg_time = np.mean(solving_times)
    median_time = np.median(solving_times)
    total_time = np.sum(solving_times)/3600

    print(colored("\nSummary:", "cyan", attrs=["bold"]))
    print(colored(f"Total files processed: {total_files}", "cyan"))
    print(colored(f"Successfully solved files: {solved_files}", "green"))
    print(colored(f"Failed files: {total_files - solved_files}", "red"))
    print(colored(f"Success percentage: {success_percentage:.2f}%", "cyan"))

    print(colored("\nSolving Time Statistics:", "yellow", attrs=["bold"]))
    print(colored(f"Maximum Time: {max_time:.2f} seconds", "yellow"))
    print(colored(f"Minimum Time: {min_time:.2f} seconds", "yellow"))
    print(colored(f"Average Time: {avg_time:.2f} seconds", "yellow"))
    print(colored(f"Median Time: {median_time:.2f} seconds", "yellow"))
    print(colored(f"Total Time: {total_time:.2f} hours", "yellow"))
else:
    print(colored("No files processed.", "red"))

print(colored("Processing complete.", "green"))
