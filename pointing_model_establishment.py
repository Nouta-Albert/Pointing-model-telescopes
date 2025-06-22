import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, AltAz, EarthLocation
from astropy.time import Time
import os
import subprocess
import glob
import sys
from termcolor import colored
from Tpoint_emulation_full import run_pointing_model

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#==========================================
#USER INPUT PARAMS - MODIFY THESE AS NEEDED
#==========================================
""" parameters to insert before running the script. Also ensure that the correct ref pixels/fiber coords are written in the header of the files, if not replace the appropriate coordinates below (x_ref and y_ref in the process_fits_images function) before running this script."""
                                    
input_dir = "/home/astro/alberte2/NewSeq-2025-06-09_non-truncated/" #directory where the solved files are located. Output files will also be saved to this diretory.
extension = "new"#replace with solved file correct ending ex .new
date = "12-06-2025" # Enter date in format DD-MM-YYY
selected_terms = ["AN", "AW", "IA", "IE", "CA", "TF", "PZZ3", "PZZ5", "HZCZ4"] #Terms to use for pointing model. For list  of available terms check Tpoint_emulation_full.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


output_dir = f"{input_dir}/Pointing_model_results"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
pointing_data_file = f"{output_dir}/PM-{date}.dat" 

# Observatory location - La Silla
obs_location = EarthLocation(lon=-70.732951945082, lat=-29.259497822801, height=2378.0)

# Step 1: Process FITS Files
def get_key(header, keys):
    for key in keys:
        value = header.get(key, None)
        if value is not None:
            return value
    return None

def process_fits_images(files_dir):
    pointing_data = []

    for file_name in os.listdir(files_dir):
        if file_name.endswith(extension): 
            fits_path = os.path.join(files_dir, file_name)

            try:
                with fits.open(fits_path) as hdul:
                    header = hdul[0].header
                    w = WCS(header)

                    # Fallback keyword options (OGE or ESO variants)
                    time_keys = [
                        'HIERARCH ESO TEL TARG TUNIX START', #Ecam new system 
                        'HIERARCH OGE TEL TARG TUNIX INST', #for ECAM old system
                        'HIERARCH ESO TEL TARG TUNIX INST'#for coralie
                    ]
                    az_keys = [
                        'HIERARCH ESO TEL TARG AZI START', #Ecam new system
                        'HIERARCH OGE TEL TARG AZI INST', #for ECAM old system
                        'HIERARCH ESO TEL TARG AZI INST' #for Coralie
                    ]     
                    el_keys = [
                        'HIERARCH ESO TEL TARG ELE START', #Ecam new system
                        'HIERARCH OGE TEL TARG ELE INST', #for Ecam old system
                        'HIERARCH ESO TEL TARG ELE INST' #for Coralie
                    ]    
                    #RA_keys = [
                    #    'RA ']
                    #DEC_keys = [
                    #    'DEC']  
                    
                    x_ref_keys = [
                        'HIERARCH OGE DET PIXELREF X',
                        'HIERARCH FIBER XREFCUR'
                    ]
                    y_ref_keys = [
                        'HIERARCH OGE DET PIXELREF Y',
                        'HIERARCH FIBER YREFCUR'
                    ]

                    # Read values using fallbacks
                    obs_time_val = get_key(header, time_keys)
                    az_expected = get_key(header, az_keys)
                    el_expected = get_key(header, el_keys)
                    #ra_expected = get_key(header, RA_keys)
                    #dec_expected = get_key(header, DEC_keys)
                    x_ref = get_key(header, x_ref_keys)
                    y_ref = get_key(header, y_ref_keys)

                    if None in [obs_time_val, az_expected, el_expected, x_ref, y_ref]:
                        print(colored(f"Missing data in {file_name}. Skipping...", "yellow"))
                        continue
                
                    obs_time = Time(float(obs_time_val), format='unix')
                    az_expected = float(az_expected) % 360
                    el_expected = float(el_expected)
                    x_ref = float(x_ref)  #2121
                    y_ref = float(y_ref)  #2121.5

                    # Convert pixel to world coords
                    pixel_world_coords = w.pixel_to_world(x_ref, y_ref)
                    sky_coord = SkyCoord(
                        ra=pixel_world_coords.ra, dec=pixel_world_coords.dec,
                        frame="icrs", unit="deg"
                    )
                    altaz_frame = AltAz(location=obs_location, obstime=obs_time)
                    altaz_coord = sky_coord.transform_to(altaz_frame)
                    
                    az_actual = (altaz_coord.az.deg + 180) % 360
                    el_actual = altaz_coord.alt.deg

                    #convert catalogue ra, dec to alt-az
                    #sky_coord2 = SkyCoord(
                     #   ra=ra_expected, dec=dec_expected,
                      #  frame="icrs", unit="deg")
                    #altaz_coord2 = sky_coord2.transform_to(altaz_frame)
                    #az_expected = (altaz_coord2.az.deg + 180) % 360
                    #el_expected = altaz_coord2.alt.deg
                    
                    # Save data
                    pointing_data.append({
                        'az_expected': az_expected,
                        'el_expected': el_expected,
                        'az_actual': az_actual,
                        'el_actual': el_actual,
                    })

            except Exception as e:
                print(colored(f"Error processing {file_name}: {e}", "red"))

    return pointing_data

    
# Step 2: Generate Data File
def create_pointing_data_file(data, pointing_file):
    try:
        with open(pointing_file, "w+") as out:
            for entry in data:
                out.write(f"{entry['az_expected']} {entry['el_expected']} {entry['az_actual']} {entry['el_actual']}\n")
        print(colored(f"Pointing data file created at {pointing_file}", "cyan"))
    except Exception as e:
        print(colored(f"Error creating pointing data file: {e}", "red"))


# Step 3: Compute Pointing Model
pointing_data = process_fits_images(input_dir)
create_pointing_data_file(pointing_data, pointing_data_file)
run_pointing_model(pointing_data_file, output_dir, selected_terms)

 