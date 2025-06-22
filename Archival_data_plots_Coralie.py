# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 11:08:31 2025

@author: Albert Einstein
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime

# Load CORALIE data from CSV files
coralie_data = pd.read_csv(r"C:\Users\Albert Einstein\Desktop\Einstein docs\Master 2 astro\Results_CORALIE_archives_analysis\pointing_data.csv")
coralie_error_data = pd.read_csv(r"C:\Users\Albert Einstein\Desktop\Einstein docs\Master 2 astro\Results_CORALIE_archives_analysis\pointing_errors.csv")

# Remove CORALIE data points with az or alt errors greater than 350 arcsecs
mask = (coralie_error_data['az_error'].abs() <= 350) & (coralie_error_data['el_error'].abs() <= 350)
coralie_data = coralie_data[mask].reset_index(drop=True)
coralie_error_data = coralie_error_data[mask].reset_index(drop=True)

# Negate the altitude errors before saving the cleaned data
coralie_error_data['el_error'] = -coralie_error_data['el_error']

# Save cleaned CORALIE data to new CSV files
coralie_data.to_csv(r"C:\Users\Albert Einstein\Desktop\Einstein docs\Master 2 astro\Results_CORALIE_archives_analysis\pointing_data_cleaned.csv", index=False)
coralie_error_data.to_csv(r"C:\Users\Albert Einstein\Desktop\Einstein docs\Master 2 astro\Results_CORALIE_archives_analysis\pointing_errors_cleaned.csv", index=False)

# Convert observation dates to datetime objects
observation_Dates = [datetime.datetime.fromisoformat(date) for date in coralie_data['Obs_date']]

ez_expected = coralie_data['az_expected'].to_numpy(dtype=float)
el_expected = coralie_data['el_expected'].to_numpy(dtype=float)
ez_actual = coralie_data['az_actual'].to_numpy(dtype=float)
el_actual = coralie_data['el_actual'].to_numpy(dtype=float)

# Extract pointing errors
az_errors = coralie_error_data['az_error'].to_numpy(dtype=float)
el_errors = coralie_error_data['el_error'].to_numpy(dtype=float)

# Plot Azimuth Error vs. Date
plt.figure(figsize=(12, 5))
plt.scatter(observation_Dates, az_errors, color='green', alpha=0.2, label="Azimuth Error [arcsec]")
plt.xlabel("Observation Date")
plt.ylabel("Azimuth Error [arcsec]")
plt.title("Pointing Error Az vs. Observation Date - CORALIE")
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=3))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
plt.xticks(rotation=45)
plt.grid()
plt.show()

# Plot Elevation Error vs. Date
plt.figure(figsize=(12, 5))
plt.scatter(observation_Dates, el_errors, color='blue', alpha=0.2, label="Elevation Error [arcsec]")
plt.xlabel("Observation Date")
plt.ylabel("Elevation Error [arcsec]")
plt.title("Pointing Error Elevation vs. Observation Date - CORALIE")
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=3))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
plt.xticks(rotation=45)
plt.grid()
plt.show()

# Define function to update colorbar with dates
def update_colorbar_with_dates(cbar, observation_dates):
    cbar_ticks = cbar.get_ticks()
    cbar_labels = [pd.to_datetime(t, unit='s').strftime('%Y-%m-%d') for t in cbar_ticks]
    cbar.set_ticklabels(cbar_labels)

# Plot Azimuth Error vs. Azimuth
plt.figure(figsize=(12, 5))
scatter = plt.scatter(ez_actual, az_errors, c=[t.timestamp() for t in observation_Dates], cmap='viridis', alpha=0.5)
plt.title('Azimuth Error vs Azimuth - CORALIE')
plt.xlabel('Azimuth [degrees]')
plt.ylabel('Azimuth Error [arcsecs]')
cbar = plt.colorbar(scatter)
update_colorbar_with_dates(cbar, observation_Dates)
cbar.set_label("Obs date")
plt.grid()
plt.show()

# Plot Azimuth Error vs Elevation
plt.figure(figsize=(12, 5))
scatter = plt.scatter(el_actual, az_errors, c=[t.timestamp() for t in observation_Dates], cmap='brg', alpha=0.5)
plt.title('Azimuth Error vs Elevation - CORALIE')
plt.xlabel('Elevation [degrees]')
plt.ylabel('Azimuth Error [arcsecs]')
cbar = plt.colorbar(scatter)
update_colorbar_with_dates(cbar, observation_Dates)
cbar.set_label("Obs date")
plt.grid()
plt.show()

# Plot Elevation Error vs Azimuth
plt.figure(figsize=(12, 5))
scatter = plt.scatter(ez_actual, el_errors, c=[t.timestamp() for t in observation_Dates], cmap='nipy_spectral', alpha=0.5)
plt.title('Elevation Error vs Azimuth - CORALIE')
plt.xlabel('Azimuth [degrees]')
plt.ylabel('Elevation Error [arcsecs]')
cbar = plt.colorbar(scatter)
update_colorbar_with_dates(cbar, observation_Dates)
cbar.set_label("Obs date")
plt.grid()
plt.show()

# Plot Elevation Error vs Elevation
plt.figure(figsize=(12, 5))
scatter = plt.scatter(el_actual, el_errors, c=[t.timestamp() for t in observation_Dates], cmap='gist_rainbow', alpha=0.5)
plt.title('Elevation Error vs Elevation - CORALIE')
plt.xlabel('Elevation [degrees]')
plt.ylabel('Elevation Error [arcsecs]')
cbar = plt.colorbar(scatter)
update_colorbar_with_dates(cbar, observation_Dates)
cbar.set_label("Obs date")
plt.grid()
plt.show()