# Pointing Model

This repository contains a comprehensive collection of Python scripts for astronomical data analysis, focusing on telescope pointing accuracy, astrometric calibration and image analysis. The tools were developed primarily for the CORALIE and ECAM instruments at La Silla Observatory but can be adapted for other telescopes.

## Script descriptions

### I. Pointing Model Development

**Tpoint_emulation_full.py**

The core script implements a complete Python-based telescope pointing calibration and modeling framework that emulates the industry-standard TPOINT system. It supports over 20 different correction terms (IA, IE, CA, TF, etc.) and performs least-squares fitting of model parameters via the TRF algorithm. The script generates comprehensive diagnostic plots including residual analysis and parameter correlation matrices, then packages all results into an interactive HTML report.

**pointing_model_establishment.py**

The complete pipeline for establishing new pointing models from raw observations. This script processes directories of solved FITS images, extracts pointing information from headers, converts coordinates and generates the input files needed for model fitting. It includes special handling for different telescope header formats and automatically integrates with Tpoint_emulation_full.py for end-to-end processing.

**Tpoint_emulation_tests.py**

This validation script compares two different pointing model configurations on the same dataset. It computes RMS errors for each model and generates comparison visualizations including polar difference plots and 2D error maps. The script outputs detailed Excel reports with numerical comparisons that are essential for model selection and tuning.

#### Usage Instructions

##### Prepare Your FITS Images
Ensure they are solved and contain the required keywords in the headers:
- Az/El coordinates
- Unix timestamp
- Detector reference pixel coordinates

##### Configure Parameters
In `pointing_model_establishment.py`, set:
- `input_dir`: directory containing solved FITS files
- `extension`: file extension (e.g., "new")
- `date`: for output naming
- `selected_terms`: terms to include in the pointing model (check `Tpoint_emulation_full.py` for full list of available terms)

##### Run the Pipeline
Execute `pointing_model_establishment.py` and inspect outputs in the generated `/Pointing_model_results` subdirectory.

#### System Components

##### 1. pointing_model_establishment.py

This is the main script that scans a directory of solved FITS images, extracts expected and actual telescope pointing (azimuth and elevation), creates a formatted data file for model fitting, and launches the pointing model computation.

#### Key Functions:

**process_fits_images()**
Processes FITS files and extracts pointing data by:
- Scanning ".new" (or user-defined) FITS files in the input directory
- Extracting expected az/el from telescope headers
- Computing actual pointing via WCS transformation of pixel coordinates (x_ref, y_ref)
- Converting actual sky positions to AltAz frame using time and observatory location
- Returning a list of structured pointing entries

**create_pointing_data_file()**
Creates input file for the pointing model by converting the processed pointing data into a plain text ".dat" file of az/el values (expected vs actual) required for model fitting.

##### 2. Tpoint_emulation_full.py

This script contains the core mathematical functions that load the pointing data file, fit a pointing model using least-squares optimization, analyze fit quality and parameter correlations, and generate visual plots with results saved in HTML output.

#### Key Functions:

**load_data(file_path)**
Loads az/el data from the ".dat" file, converting degrees to radians for internal calculations.

**pointing_model(params, az_real, el_real, selected_terms)**
Computes pointing corrections (ΔAz, ΔEl) using user-selected model terms (e.g., AN, AW, IA, etc.) based on standard astronomical pointing model formulations.

**residuals(params, az_real, el_real, az_actual, el_actual, selected_terms)**
Calculates the residuals between observed and model-corrected az/el coordinates for optimization algorithms.

**fit_pointing_model(data, selected_terms, initial_params)**
Performs least-squares fitting using `scipy.optimize.least_squares` and returns:
- Best-fit parameter values
- Parameter uncertainties (converted to arcseconds)
- Covariance matrix for correlation analysis

**calculate_rms(data, fitted_params, selected_terms)**
Computes RMS residuals in azimuth and elevation in arcseconds to quantify model performance.

**analyze_correlations(cov_matrix, selected_terms)**
Builds a correlation matrix from the covariance and identifies parameter pairs with correlation |r| > 0.8 to detect potential degeneracies.

**analyze_and_plot_all(...)**
Creates a comprehensive 3×3 diagnostic grid plot including:
- Sky distribution of data points
- Residuals in azimuth and elevation
- Correlation heatmaps
- Histograms of residuals and scatter diagnostics

**save_results(...)**
Saves fitted parameter values with uncertainties, correlation matrix, and diagnostic plot grid, then generates a comprehensive HTML report.

**run_pointing_model(...)**
The public interface for external scripts that triggers the full analysis and saving workflow.

#### Results - Output Files

The system produces an HTML file in the specified output directory containing:

##### Numerical Results
- Fitted parameter values with 1-sigma uncertainties
- RMS fit errors for azimuth and elevation
- Information about strongly correlated parameters

##### Correlation Matrix
Complete matrix showing correlations between all parameters, helpful for detecting degeneracy or overfitting issues.

##### Visualizations
A 3×3 grid of diagnostic plots showing:
- Point distribution on the sky
- Residuals of the fitting process in Az and El
- Parameter correlation heatmap
- Histograms of residuals for quality assessment

-----------------------------------------------------------------------------------------------------------------------

### II. Data processing and analysis
**Coralie_analysis_script.py / Ecam_analysis_script.py**  
These automated pipelines process telescope observation files to calculate pointing errors. They handle the complete workflow from plate solving FITS images to computing expected versus actual pointing positions. The scripts calculate azimuth and elevation errors while handling instrument-specific header formats and log solving statistics to CSV files. The ECAM version includes more robust error handling for diverse observation types.

**Archival data plots Coralie.py / Archival data plots Ecam.py**  
These scripts visualize historical pointing error data from telescope archives. They load and clean pointing data from CSV files, generate time-series plots of azimuth/elevation errors and create scatter plots of errors versus telescope position. The ECAM version includes additional weekly error statistics and more sophisticated outlier removal.

**Pointing error coralie analysis.ipynb**  
This analysis notebook focuses specifically on quantifying CORALIE's pointing errors. It loads pre-processed CSV data, converts pixel offsets to sky offsets and calculates pointing errors in both azimuth and elevation. The notebook generates multiple visualization formats including time-series plots, polar plots, and error distribution histograms to help identify systematic pointing issues.

**Pointing error Coralie analysis using plate solving.ipynb**  
Similar to the ECAM version but tailored for CORALIE observations, this notebook handles the plate solving workflow for CORALIE's specific FITS format. It includes additional header parsing logic for CORALIE's metadata structure and outputs detailed logs of the solving process with timestamps for performance monitoring.

**Pointing error ECAM analysis using my plate solving.ipynb**  
This Jupyter Notebook performs plate solving on ECAM FITS images in an interactive environment. It systematically processes each FITS file by extracting RA and Dec from headers, executes the solve-field command and logs the results. The notebook format allows for real-time inspection of solving results and immediate troubleshooting when plate solving fails.

**stars_elongation.py**  
A specialized tool for analyzing stellar elongation patterns that may indicate optical or tracking issues. The script detects stars, measures their PSF shapes and orientations and overlays the expected alt-az grid based on derotator angle. It produces multi-panel visualizations showing both the full image context and detailed star analyses.

-----------------------------------------------------------------------------------------------------------------------

### III. Instrument calibration and verification

**Compute_transformation_matrix_coralie.ipynb**  
This interactive notebook calculates the critical transformation between CORALIE's detector coordinates and telescope pointing coordinates. It measures star positions with known telescope offsets and computes the transformation matrix. The notebook format allows step-by-step validation of the transformation calculations.

**Solvefield_accuracy_test.py**  
This verification script rigorously tests astrometric solution accuracy by cross-matching detected stars with Gaia DR2 catalog positions. It performs comprehensive error analysis including separation histograms and sector-based error mapping to identify systematic trends. The script generates multiple visualization formats to help diagnose plate solving quality across the entire field of view.

**Rotation_matrix_criteria.py**  
A technical validation tool for WCS solutions in FITS files. The script analyzes CD matrices through polar decomposition to separate rotation from scaling components. 

