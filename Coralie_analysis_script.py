import os
import subprocess
import pandas as pd
from datetime import datetime
from astropy.io import fits
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, ICRS
from astropy.time import Time
from astropy.wcs import WCS
import numpy as np
import shutil
import glob
import time
from termcolor import colored
from datetime import timedelta
from astropy import units as u

##############################################################################
# Directories
archives_dir = "/data/CORALIE14RAW/raw/"
temp_working_dir = "/home/astro/alberte2/Archives_CORALIE/"
solved_output_dir = os.path.join(temp_working_dir, "solved/")
os.makedirs(solved_output_dir, exist_ok=True)
output_dir = "/home/astro/alberte2/Results_CORALIE_archives_analysis/"
##############################################################################

# Output files
unsolved_files_log = os.path.join(output_dir, "unsolved_files.txt")
solving_stats_file = os.path.join(output_dir, "solving_stats.csv")
pointing_errors_file = os.path.join(output_dir, "pointing_errors.csv")
pointing_data_file = os.path.join(output_dir, "pointing_data.csv")

# Observatory location
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

# Helper function to get fiber coordinates
def get_fiber_coordinates(file_date):
    if datetime(2014, 11, 13) < file_date < datetime(2021, 12, 14):
        return 652, 380
    elif datetime(2021, 12, 14) <= file_date < datetime(2024, 1, 1):
        return 653, 384
    elif file_date >= datetime(2024, 1, 1):
        return 652, 380
    else:
        return 652, 380
        
def get_science_files(start_date, end_date):
    # Convert user-defined dates to datetime objects
    start_dt = datetime.strptime(start_date, '%Y-%m-%d')
    end_dt = datetime.strptime(end_date, '%Y-%m-%d')

    fits_files = []
    for root, _, files in os.walk(archives_dir):
        for file in files:
            if "PCORALIE" in file and file.endswith(".fits"):
                try:
                    # Extract observation date from the filename (e.g., "2025-02-07T03:39:42")
                    obs_date_str = file.split('.')[1]  # Extract the "2025-02-07T03:39:42" part
                    obs_date = datetime.strptime(obs_date_str.split('T')[0], '%Y-%m-%d')  # Only "YYYY-MM-DD"

                    # Check if the observation date falls within the user-defined range
                    if start_dt <= obs_date <= end_dt:
                        fits_files.append(os.path.join(root, file))

                except ValueError:
                    print(f"Skipping file with invalid date format: {file}")

    # Sort the files by observation date in reverse order (most recent first)
    fits_files.sort(
        key=lambda x: datetime.strptime(x.split('.')[1], '%Y-%m-%dT%H:%M:%S'),
        reverse=True
    )
    return fits_files    

# Function to perform platesolving and log results
def plate_solve_and_log(fits_file):
    file_name = os.path.basename(fits_file)
    output_name = file_name.replace(".fits", ".solved.fits")
    solved_file_path = os.path.join(solved_output_dir, output_name)

    with fits.open(fits_file) as hdul:
        header = hdul[0].header

        # Try to get RA and DEC from header
        ra_approx = header.get('RA', None) or header.get('HIERARCH ESO TEL TARG ALPHA INST', None)
        dec_approx = header.get('DEC', None) or header.get('HIERARCH ESO TEL TARG DELTA INST', None)

        # If RA/DEC are missing, check for AZI/ELE and convert
        if ra_approx is None or dec_approx is None:
            azi = header.get('HIERARCH AZI', None)
            ele = header.get('HIERARCH ELE', None)
            unix_time = header.get('HIERARCH UNIXTIME', None)

            if azi is not None and ele is not None and unix_time is not None:
                try:
                    # Convert to floats
                    azi = float(azi)
                    ele = float(ele)
                    unix_time = float(unix_time)

                    # Convert observation time to Astropy Time object
                    obs_time = Time(unix_time, format='unix')

                    # Convert Az/El to RA/Dec
                    altaz_coord = SkyCoord(az=(azi - 180 )* u.deg, alt=ele * u.deg, frame=AltAz(location=obs_location, obstime=obs_time))
                    icrs_coord = altaz_coord.transform_to(ICRS)

                    ra_approx = icrs_coord.ra.deg
                    dec_approx = icrs_coord.dec.deg

                    print(colored(f"Converted AZI/ELE to RA/DEC: ({ra_approx:.6f}, {dec_approx:.6f})", "yellow"))

                except Exception as e:
                    print(colored(f"Error converting AZI/ELE to RA/DEC: {e}", "red"))
                    ra_approx, dec_approx = None, None  # Fallback to no RA/DEC

    # Build the `solve-field` command dynamically
    command = f"solve-field --overwrite --no-verify --scale-units arcsecperpix --scale-low 0.2 --scale-high 0.4 " \
              f"--radius 5 --depth 100 --objs 200 --no-plots --dir {solved_output_dir} --new-fits {output_name}"

    # If RA/DEC are available, add them
    if ra_approx is not None and dec_approx is not None:
        command += f" --ra {ra_approx} --dec {dec_approx}"

    # Append the FITS file at the end (IMPORTANT: `solve-field` expects it last)
    command += f" {fits_file}"

    try:
        start_time = time.time()

        # Run command in shell mode
        result = subprocess.run(command, shell=True, cwd=solved_output_dir, capture_output=True, text=True)

        # Debugging output (optional)
        # print(f"STDOUT: {result.stdout}")
        # print(f"STDERR: {result.stderr}")

        end_time = time.time()
        elapsed_time = end_time - start_time

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

# Function to find closest CORALIE file for a given PCORALIE file
def find_closest_coralie_file(pcor_path):
    pcor_dir = os.path.dirname(pcor_path)
    pcor_name = os.path.basename(pcor_path)
    pcor_date_str = pcor_name.split('.')[1]  # Extract "YYYY-MM-DDTHH:MM:SS"
    pcor_date = datetime.strptime(pcor_date_str, '%Y-%m-%dT%H:%M:%S')

    # First, try finding a CORALIE file with the exact same timestamp
    exact_match = os.path.join(pcor_dir, pcor_name.replace("PCORALIE", "CORALIE"))
    if os.path.exists(exact_match):
        print(colored(f"Using exact CORALIE match: {exact_match}", "cyan"))
        return exact_match

    # If no exact match, find the closest CORALIE file by timestamp
    coralie_files = [f for f in os.listdir(pcor_dir) if f.startswith("CORALIE") and f.endswith(".fits")]
    closest_file = None
    min_time_diff = float('inf')

    for cor_file in coralie_files:
        cor_date_str = cor_file.split('.')[1]  # Extract "YYYY-MM-DDTHH:MM:SS"
        try:
            cor_date = datetime.strptime(cor_date_str, '%Y-%m-%dT%H:%M:%S')
            time_diff = abs((pcor_date - cor_date).total_seconds())

            if time_diff < min_time_diff:
                min_time_diff = time_diff
                closest_file = os.path.join(pcor_dir, cor_file)
        except ValueError:
            continue

    if closest_file:
        print(colored(f"Closest CORALIE file found: {closest_file}", "cyan"))
    else:
        print(colored(f"No matching CORALIE file found for {pcor_path}", "red"))

    return closest_file

# Function to extract target coordinates
def get_target_coordinates(header, file_date, file_path):
    threshold_date = datetime(2023, 12, 11)
    
    if file_date < threshold_date:
        ra_str = header.get('HIERARCH ESO OBS ALPHACAT', None)
        dec_str = header.get('HIERARCH ESO OBS DELTACAT', None)
    else:
        coralie_file_path = find_closest_coralie_file(file_path)
        if coralie_file_path is None:
            print(colored(f"Missing CORALIE file for {file_path}. Skipping...", "red"))
            return None, None

        try:
            with fits.open(coralie_file_path) as coralie_hdul:
                coralie_header = coralie_hdul[0].header
                ra_str = coralie_header.get('HIERARCH ESO OBS ALPHACAT', None)
                dec_str = coralie_header.get('HIERARCH ESO OBS DELTACAT', None)
        except FileNotFoundError:
            print(colored(f"Could not open CORALIE file for {file_path}.", "red"))
            return None, None

    if ra_str and dec_str:
        coords = SkyCoord(ra_str, dec_str, unit=("hourangle", "deg"))
        return coords.ra.deg, coords.dec.deg
    return None, None

# Function to calculate pointing errors
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

                    # Extract observation time
                    obs_date_str = header.get('DATE-OBS', None)
                    if obs_date_str:
                        unix_time = header.get('HIERARCH ESO TEL TARG TUNIX INST', None)
                        if unix_time:
                            obs_time = Time(float(unix_time), format='unix')  # Use the UNIX time from the header
                        else:
                            obs_time = Time(obs_date_str).unix  # Convert DATE-OBS to UNIX time
                            obs_time = Time(obs_time, format='unix')  # Ensure proper formatting
                    else:
                        unix_time = header.get('STARTTIM', None)
                        if unix_time:
                            obs_time = Time(float(unix_time), format='unix')
                            obs_date_str = obs_time.isot
                        else:
                            print(colored(f"Missing observation time for {file_name}. Skipping...", "red"))
                            continue

                    # Extract file date from filename
                    file_date_str = file_name.split('.')[1]  # Extract "YYYY-MM-DDTHH:MM:SS"
                    file_date = datetime.strptime(file_date_str.split('T')[0], '%Y-%m-%d')

                    # Determine the original archive file path
                    original_date = file_name.split('.')[1][:10]  # Extract "YYYY-MM-DD"
                    original_file_name = file_name.replace(".solved.fits", ".fits")

                    # First, assume the file is stored in the directory matching its date
                    original_file_path = os.path.join(archives_dir, original_date, original_file_name)

                    # If file is not found, check the previous day's archive folder
                    if not os.path.exists(original_file_path):
                        previous_date = (datetime.strptime(original_date, "%Y-%m-%d") - timedelta(days=1)).strftime("%Y-%m-%d")
                        original_file_path = os.path.join(archives_dir, previous_date, original_file_name)

                    # Ensure the file exists before proceeding
                    if not os.path.exists(original_file_path):
                        print(colored(f"Original file not found for {file_name}. Skipping...", "red"))
                        continue

                    # Get expected RA/DEC coordinates from the correct file
                    ra, dec = get_target_coordinates(header, file_date, original_file_path)
                    if ra is None or dec is None:
                        continue

                    # Convert to horizontal coordinates
                    sky_coord = SkyCoord(ra, dec, unit="deg", frame="icrs")
                    altaz_frame = AltAz(location=obs_location, obstime=obs_time)
                    altaz_coord = sky_coord.transform_to(altaz_frame)

                    az_expected = altaz_coord.az.deg
                    el_expected = altaz_coord.alt.deg

                    # Retrieve fiber coordinates
                    fiber_x = header.get('HIERARCH FIBER XREFCUR', None)
                    fiber_y = header.get('HIERARCH FIBER YREFCUR', None)
                    if fiber_x is None or fiber_y is None:
                        fiber_x, fiber_y = get_fiber_coordinates(file_date)

                    # Convert fiber coordinates to horizontal coordinates
                    pixel_world_coords = w.pixel_to_world(fiber_x, fiber_y)
                    fiber_sky_coord = SkyCoord(
                        ra=pixel_world_coords.ra, dec=pixel_world_coords.dec, frame="icrs", unit="deg"
                    )
                    fiber_altaz_coord = fiber_sky_coord.transform_to(altaz_frame)

                    az_actual = fiber_altaz_coord.az.deg
                    el_actual = fiber_altaz_coord.alt.deg

                    # Calculate pointing errors
                    el_err = (el_actual - el_expected) * 3600
                    az_err = (az_actual - az_expected) * np.cos(np.radians(el_expected)) * 3600

                    # Save pointing data
                    pointing_data.append({
                        'Obs_date': obs_date_str,
                        'az_expected': az_expected,
                        'el_expected': el_expected,
                        'az_actual': az_actual,
                        'el_actual': el_actual
                    })
                    pointing_errors.append({
                        'Timestamp': obs_date_str,
                        'az_error': az_err,
                        'el_error': el_err
                    })

            except Exception as e:
                print(colored(f"Error processing {file_name}: {e}", "red"))

    # Save results
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

####################################################################################################
# User-defined start and end dates
start_date = '2014-11-13' #input(colored("Enter start date (YYYY-MM-DD): ", "cyan", attrs=["bold"]))
end_date = '2025-02-17' #input(colored("Enter end date (YYYY-MM-DD): ", "magenta", attrs=["bold"]))
####################################################################################################

# Process each FITS file
science_files = get_science_files(start_date, end_date)
for fits_file in science_files:
    total_files += 1
    print(colored(f"Platesolving: {fits_file}", "magenta"))
    temp_file_path = os.path.join(temp_working_dir, os.path.basename(fits_file))
    shutil.copy(fits_file, temp_file_path)
    
    if plate_solve_and_log(temp_file_path):
        solved_files += 1
        
        #print("Files to process for pointing errors:")
        #print(os.listdir(solved_output_dir))
        
        calculate_and_log_pointing_errors()

    solving_times.append(time.time() - os.path.getmtime(temp_file_path))
    for file in os.listdir(solved_output_dir):
        os.remove(os.path.join(solved_output_dir, file))
    
    for item in os.listdir(temp_working_dir):
        item_path = os.path.join(temp_working_dir, item)
        if item != "solved":
            if os.path.isfile(item_path):
                os.remove(item_path)
            elif os.path.isdir(item_path):
                shutil.rmtree(item_path)
                
# Calculate and display success statistics
# Read the number of entries in pointing_errors.csv
try:
    pointing_errors_df = pd.read_csv(pointing_errors_file)
    num_pointing_errors = len(pointing_errors_df)
except FileNotFoundError:
    num_pointing_errors = 0

# Compute files that solved but didn't produce pointing errors
unsolved_pointing_errors = solved_files - num_pointing_errors
unsolved_pointing_percentage = (unsolved_pointing_errors / solved_files) * 100 if solved_files > 0 else 0

# Compute statistics if there were files processed
if total_files > 0:
    success_percentage = (solved_files / total_files) * 100
    max_time = np.max(solving_times)
    min_time = np.min(solving_times)
    avg_time = np.mean(solving_times)
    total_time = np.sum(solving_times) / 3600

    print(colored("\nSummary:", "cyan", attrs=["bold"]))
    print(colored(f"Total files processed: {total_files}", "cyan"))
    print(colored(f"Successfully solved files: {solved_files}", "green"))
    print(colored(f"Failed files: {total_files - solved_files}", "red"))
    print(colored(f"Success percentage: {success_percentage:.2f}%", "magenta"))
    print(colored(f"Total Time: {total_time:.2f} hours", "yellow"))
    
    # New statistics
    print(colored(f"Solved files without pointing errors: {unsolved_pointing_errors}", "red"))
    print(colored(f"Percentage of solved files without pointing errors: {unsolved_pointing_percentage:.2f}%", "red"))

else:
    print(colored("No files processed.", "red"))

print(colored("Processing complete.", "green"))