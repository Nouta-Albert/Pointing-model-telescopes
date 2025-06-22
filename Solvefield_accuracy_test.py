import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from astropy.table import Table
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.visualization import simple_norm
import random
from astropy.visualization import ZScaleInterval, ImageNormalize

# ----------------------
# Configurable Parameters
# ----------------------
FITS_FILE = '/home/astro/alberte2/PM_test_LPM/Cross_matching/ECAM.2025-04-02T08:39:44.000.solved.fits'      
OUTPUT_PNG = '/home/astro/alberte2/PM_test_LPM/Cross_matching/detected_stars7.png'
HIST = "/home/astro/alberte2/PM_test_LPM/Cross_matching/separation_histogram7.png"
FWHM = 10                         # Typical PSF FWHM in pixels
THRESHOLD_SIGMA = 10.0             # Detection threshold in sigma
RADIUS_MATCH_ARCSEC = 2.0         # Matching radius for crossmatch (arcsec)
GAIA_CATALOG = 'I/345/gaia2'      # Gaia DR2 catalog ID in Vizier
SEARCH_RADIUS_DEG = 0.25           # 16 arcmin search radius
NUM_RANDOM_DISPLAY = 5            # Number of random matches to display

# ----------------------
# Read FITS and WCS
# ----------------------
hdul = fits.open(FITS_FILE)
data = hdul[0].data
header = hdul[0].header
wcs = WCS(header)

# ----------------------
# Detect stars
# ----------------------
mean, median, std = sigma_clipped_stats(data, sigma=3.0)
daofind = DAOStarFinder(fwhm=FWHM, threshold=THRESHOLD_SIGMA * std)
sources = daofind(data - median)

if sources is None or len(sources) == 0:
    raise RuntimeError("No stars detected.")

print(f" Detected {len(sources)} stars")

# ----------------------
# Save image with detections
# ----------------------
plt.figure(figsize=(10, 10))
zscale = ZScaleInterval()
vmin, vmax = zscale.get_limits(data)
norm = ImageNormalize(vmin=vmin, vmax=vmax)
plt.imshow(data, cmap='gray', origin='lower', norm=norm)
plt.scatter(sources['xcentroid'], sources['ycentroid'], s=30, edgecolor='red', facecolor='none')
plt.title(f"Detected Stars: {len(sources)}")
plt.savefig(OUTPUT_PNG)
plt.close()

# ----------------------
# Convert to sky coordinates using WCS
# ----------------------
sky_coords = wcs.pixel_to_world(sources['xcentroid'], sources['ycentroid'])
source_skycoord = SkyCoord(ra=sky_coords.ra, dec=sky_coords.dec)

# ----------------------
# Query Gaia catalog using header RA/DEC
# ----------------------
center = SkyCoord(ra=header['RA'] * u.deg, dec=header['DEC'] * u.deg)
search_radius = SEARCH_RADIUS_DEG * u.deg

Vizier.ROW_LIMIT = -1
result = Vizier.query_region(center, radius=search_radius, catalog=GAIA_CATALOG)

if not result:
    raise ValueError("No Gaia sources found in the field of view.")

gaia_table = result[0]
gaia_coords = SkyCoord(ra=gaia_table['RA_ICRS'], dec=gaia_table['DE_ICRS'])

print(f" Retrieved {len(gaia_coords)} Gaia sources in field")

# ----------------------
# Crossmatch
# ----------------------
idx, d2d, _ = source_skycoord.match_to_catalog_sky(gaia_coords)
matched = d2d.arcsec < RADIUS_MATCH_ARCSEC
n_matched = np.sum(matched)

print(f" Crossmatched {n_matched} out of {len(source_skycoord)} detected sources")

# ----------------------
# Display random matched pairs
# ----------------------
matched_indices = np.where(matched)[0]

print("\n Sample matched sources and separations (arcsec):")
for i in random.sample(list(matched_indices), min(NUM_RANDOM_DISPLAY, n_matched)):
    print(f"Detected RA/Dec: {source_skycoord[i].ra.deg:.6f}, {source_skycoord[i].dec.deg:.6f} | "
          f"Gaia RA/Dec: {gaia_coords[idx[i]].ra.deg:.6f}, {gaia_coords[idx[i]].dec.deg:.6f} | "
          f"Separation: {d2d[i].arcsec:.3f}\"")

# ----------------------
# Plot histogram of separations
# ----------------------
plt.figure(figsize=(8, 5))
plt.hist(d2d[matched].arcsec, bins=30, color='steelblue', histtype = 'step')
plt.xlabel('Separation (arcsec)')
plt.ylabel('Number of Matches')
plt.title('Histogram of Crossmatch Separations')
plt.tight_layout()
plt.savefig(HIST)
plt.show()

# ----------------------
# Plot separation vs distance from image center
# ----------------------
# Image center (assuming square-ish image)
x_center = 2121
y_center = 2121.5

# Only keep matched stars
x_matched = sources['xcentroid'][matched]
y_matched = sources['ycentroid'][matched]
sep_matched = d2d[matched].arcsec

# Compute distance from image center in pixels
radii_pixels = np.sqrt((x_matched - x_center)**2 + (y_matched - y_center)**2)




# Angular sector analysis
import matplotlib.gridspec as gridspec

# Plot separation vs radius in subplots (one per angular sector)
# Angular sector analysis
# Angular sector analysis
n_sectors = 16  # divide 360° into 8 sectors (every 45°)
sector_width_deg = 360 / n_sectors

# Compute angles from image center (in degrees, range 0-360)
angles_deg = (np.degrees(np.arctan2(y_matched - y_center, x_matched - x_center)) + 360) % 360
sector_indices = np.floor(angles_deg / sector_width_deg).astype(int)

n_cols = 4
n_rows = int(np.ceil(n_sectors / n_cols))
fig = plt.figure(figsize=(4 * n_cols, 3.5 * n_rows))
gs = gridspec.GridSpec(n_rows, n_cols)
colors = plt.cm.tab10(np.linspace(0, 1, n_sectors))

for i in range(n_sectors):
    ax = fig.add_subplot(gs[i])
    in_sector = sector_indices == i
    if np.sum(in_sector) == 0:
        continue
    ax.scatter(radii_pixels[in_sector], sep_matched[in_sector],
               s=15, alpha=0.5, color=colors[i])
    angle_min = int(i * sector_width_deg)
    angle_max = int((i + 1) * sector_width_deg)
    ax.set_title(f"Sector {angle_min}°–{angle_max}°", fontsize=10)
    ax.set_xlabel("Radius (px)", fontsize=9)
    ax.set_ylabel("Sep. (arcsec)", fontsize=9)
    ax.tick_params(axis='both', labelsize=8)
    ax.grid(True)

plt.tight_layout()
plt.suptitle("Separation vs Distance by Angular Sector", fontsize=14, y=1.02)
plt.subplots_adjust(top=0.92)
plt.savefig("/home/astro/alberte2/PM_test_LPM/Cross_matching/separation_by_sector_subplots7.png")
plt.show()


