# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 15:58:11 2025

@author: Albert Einstein
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime

# Load data from CSV files
ecam_data = pd.read_csv(r"C:\Users\Albert Einstein\Desktop\Einstein docs\Master 2 astro\Results_ECAM_archives_analysis\pointing_data.csv")
error_data = pd.read_csv(r"C:\Users\Albert Einstein\Desktop\Einstein docs\Master 2 astro\Results_ECAM_archives_analysis\pointing_errors.csv")

# Function to remove outliers while maintaining correspondence
def remove_outliers(df1, df2, column):
    Q1 = df2[column].quantile(0.25)
    Q3 = df2[column].quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 10 * IQR
    upper_bound = Q3 + 10 * IQR
    mask = (df2[column] >= lower_bound) & (df2[column] <= upper_bound)
    return df1[mask].reset_index(drop=True), df2[mask].reset_index(drop=True)

# Removing outliers while maintaining order alignment
ecam_data, error_data = remove_outliers(ecam_data, error_data, 'az_error')
ecam_data, error_data = remove_outliers(ecam_data, error_data, 'el_error')

# Convert observation dates to datetime objects
observation_Dates = [datetime.datetime.fromisoformat(date) for date in ecam_data['Obs_date']]

ez_expected = ecam_data['az_expected'].to_numpy(dtype=float)
el_expected = ecam_data['el_expected'].to_numpy(dtype=float)
ez_actual = ecam_data['az_actual'].to_numpy(dtype=float)
el_actual = ecam_data['el_actual'].to_numpy(dtype=float)

# Extract pointing errors
az_errors = error_data['az_error'].to_numpy(dtype=float)
el_errors = error_data['el_error'].to_numpy(dtype=float)

# Plot Azimuth Error vs. Date
plt.figure(figsize=(12, 5))
plt.scatter(observation_Dates, az_errors, color='green', alpha=0.1, label="Azimuth Error [arcsec]")
plt.xlabel("Observation Date")
plt.ylabel("Azimuth Error [arcsec]")
plt.title("Pointing Error Az vs. Observation Date")
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=5))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
plt.xticks(rotation=45)
plt.grid()
plt.show()

# Plot Elevation Error vs. Date
plt.figure(figsize=(12, 5))
plt.scatter(observation_Dates, el_errors, color='blue', alpha=0.1, label="Elevation Error [arcsec]")
plt.xlabel("Observation Date")
plt.ylabel("Elevation Error [arcsec]")
plt.title("Pointing Error Elevation vs. Observation Date")
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=5))
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
scatter = plt.scatter(ez_actual, az_errors, c=[t.timestamp() for t in observation_Dates], cmap='viridis', alpha=0.2)
plt.title('Azimuth Error vs Azimuth')
plt.xlabel('Azimuth [degrees]')
plt.ylabel('Azimuth Error [arcsecs]')
cbar = plt.colorbar(scatter)
update_colorbar_with_dates(cbar, observation_Dates)
cbar.set_label("Obs date")
plt.grid()
plt.show()

# Plot Azimuth Error vs Elevation
plt.figure(figsize=(12, 5))
scatter = plt.scatter(el_actual, az_errors, c=[t.timestamp() for t in observation_Dates], cmap='brg', alpha=0.2)
plt.title('Azimuth Error vs Elevation')
plt.xlabel('Elevation [degrees]')
plt.ylabel('Azimuth Error [arcsecs]')
cbar = plt.colorbar(scatter)
update_colorbar_with_dates(cbar, observation_Dates)
cbar.set_label("Obs date")
plt.grid()
plt.show()

# Plot Elevation Error vs Azimuth
plt.figure(figsize=(12, 5))
scatter = plt.scatter(ez_actual, el_errors, c=[t.timestamp() for t in observation_Dates], cmap='nipy_spectral', alpha=0.2)
plt.title('Elevation Error vs Azimuth')
plt.xlabel('Azimuth [degrees]')
plt.ylabel('Elevation Error [arcsecs]')
cbar = plt.colorbar(scatter)
update_colorbar_with_dates(cbar, observation_Dates)
cbar.set_label("Obs date")
plt.grid()
plt.show()

# Plot Elevation Error vs Elevation
plt.figure(figsize=(12, 5))
scatter = plt.scatter(el_actual, el_errors, c=[t.timestamp() for t in observation_Dates], cmap='gist_rainbow', alpha=0.2)
plt.title('Elevation Error vs Elevation')
plt.xlabel('Elevation [degrees]')
plt.ylabel('Elevation Error [arcsecs]')
cbar = plt.colorbar(scatter)
update_colorbar_with_dates(cbar, observation_Dates)
cbar.set_label("Obs date")
plt.grid()
plt.show()





import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime

# Load data from CSV files
ecam_data = pd.read_csv(r"C:\Users\Albert Einstein\Desktop\Einstein docs\Master 2 astro\Results_ECAM_archives_analysis\pointing_data.csv")
error_data = pd.read_csv(r"C:\Users\Albert Einstein\Desktop\Einstein docs\Master 2 astro\Results_ECAM_archives_analysis\pointing_errors.csv")

# Function to remove outliers while maintaining correspondence
def remove_outliers(df1, df2, column):
    Q1 = df2[column].quantile(0.25)
    Q3 = df2[column].quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 10 * IQR
    upper_bound = Q3 + 10 * IQR
    mask = (df2[column] >= lower_bound) & (df2[column] <= upper_bound)
    return df1[mask].reset_index(drop=True), df2[mask].reset_index(drop=True)

# Removing outliers while maintaining order alignment
ecam_data, error_data = remove_outliers(ecam_data, error_data, 'az_error')
ecam_data, error_data = remove_outliers(ecam_data, error_data, 'el_error')

# Convert observation dates to datetime objects
observation_Dates = [datetime.datetime.fromisoformat(date) for date in ecam_data['Obs_date']]

# Extract pointing errors
az_errors = error_data['az_error'].to_numpy(dtype=float)
el_errors = error_data['el_error'].to_numpy(dtype=float)

# Convert to DataFrame for easy resampling
df = pd.DataFrame({
    'date': observation_Dates,
    'az_error': az_errors,
    'el_error': el_errors
}).set_index('date')

# Define robust aggregation functions that handle empty groups
def safe_percentile(x, p):
    if len(x) == 0:
        return np.nan
    return np.nanpercentile(x, p)

def safe_mean(x):
    if len(x) == 0:
        return np.nan
    return np.nanmean(x)

# Resample by week and compute statistics
df_weekly = df.resample('W').agg({
    'az_error': [
        ('mean', safe_mean),
        ('p10', lambda x: safe_percentile(x, 10)),
        ('p90', lambda x: safe_percentile(x, 90))
    ],
    'el_error': [
        ('mean', safe_mean),
        ('p10', lambda x: safe_percentile(x, 10)),
        ('p90', lambda x: safe_percentile(x, 90))
    ]
})

# Flatten column names
df_weekly.columns = ['_'.join(col).strip() for col in df_weekly.columns.values]

# Drop weeks with no data
df_weekly = df_weekly.dropna()

# Plot Azimuth
plt.figure(figsize=(14, 6))
plt.fill_between(df_weekly.index, 
                df_weekly['az_error_p10'], 
                df_weekly['az_error_p90'], 
                color='skyblue', alpha=0.3, label='10th-90th Percentile')
plt.plot(df_weekly.index, df_weekly['az_error_mean'], 
         color='navy', label='Weekly Mean', lw=2)
plt.title("Azimuth Error Variation Over Time", fontsize=14)
plt.xlabel("Date", fontsize=12)
plt.ylabel("Error [arcsec]", fontsize=12)
plt.legend(fontsize=12)
plt.grid(alpha=0.3)

# Format x-axis
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=3))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# Plot Elevation
plt.figure(figsize=(14, 6))
plt.fill_between(df_weekly.index, 
                df_weekly['el_error_p10'], 
                df_weekly['el_error_p90'], 
                color='lightcoral', alpha=0.3, label='10th-90th Percentile')
plt.plot(df_weekly.index, df_weekly['el_error_mean'], 
         color='darkred', label='Weekly Mean', lw=2)
plt.title("Elevation Error Variation Over Time", fontsize=14)
plt.xlabel("Date", fontsize=12)
plt.ylabel("Error [arcsec]", fontsize=12)
plt.legend(fontsize=12)
plt.grid(alpha=0.3)

# Format x-axis
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=3))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()