import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import ZScaleInterval
import astropy.units as u
from photutils.detection import DAOStarFinder
from photutils.aperture import EllipticalAperture
from photutils.morphology import data_properties
from astropy.stats import sigma_clipped_stats


###########################################################################################################################
# Input parameters-adjust at your convenience
input_fits = "/home/astro/alberte2/euler02_data_ECAM/2025-01-31/ECAM.2025-02-01T07:37:31.000.fits"
output_dir = "/home/astro/alberte2/Ecam_alt_az_grid/"
fwhm = 13  # FWHM for star detection
threshold = 14  # Threshold for star detection in units of background std
random_seed = None  # Seed for reproducibility
r = 1.7 #ellipse scaling ratio
size = 20 #crop size
############################################################################################################################


num_stars_to_analyze = 12  
def plot_with_grids_and_stars(solved_fits, output_dir, num_stars_to_analyze, fwhm, threshold, random_seed):
    # Open the FITS file
    with fits.open(solved_fits) as hdul:
        header = hdul[0].header
        wcs = WCS(header)
        image_data = hdul[0].data

        # Extract reference pixel (optical axis)
        ref_x = header.get("HIERARCH OGE DET PIXELREF X")
        ref_y = header.get("HIERARCH OGE DET PIXELREF Y")

        # Extract rotation angle (derotator angle)
        rot_angle = -header.get("HIERARCH OGE TEL TARG ROT INST") * u.deg

        # Compute image dimensions
        ny, nx = image_data.shape

        # Display image with zscale
        zscale = ZScaleInterval()
        vmin, vmax = zscale.get_limits(image_data)
        fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': wcs})
        ax.imshow(image_data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)

        # Plot celestial grid with 10 grid lines per axis in a striking color
        ax.coords.grid(True, color='yellow', linestyle='dashed', alpha=0.7)
        ax.coords[0].set_axislabel("Right Ascension", color='black', fontsize=14)
        ax.coords[1].set_axislabel("Declination", color='black', fontsize=14)
        ax.coords[0].set_ticks(number=6)
        ax.coords[1].set_ticks(number=6)

        # Define slopes for altitude and azimuth lines
        m_alt = np.tan(rot_angle.to(u.rad).value)
        m_az = -1 / m_alt  # Perpendicular slope

        # Compute intercepts (c = y - m*x)
        c_alt = ref_y - m_alt * ref_x
        c_az = ref_y - m_az * ref_x

        # Function to find the valid segment of the line inside the image
        def get_line_points(m, c, nx, ny):
            """ Finds the valid x, y points for a line within image bounds """
            points = []

            # Intersect with left and right borders (x = 0 and x = nx)
            for x in [0, nx]:
                y = m * x + c
                if 0 <= y <= nx:
                    points.append((x, y))

            # Intersect with top and bottom borders (y = 0 and y = ny)
            for y in [0, ny]:
                if m != 0:  # Avoid division by zero
                    x = (y - c) / m
                    if 0 <= x <= ny:
                        points.append((x, y))

            return np.array(points).T  # Transpose to separate x and y

        # Get valid segments for altitude and azimuth lines
        alt_x, alt_y = get_line_points(m_alt, c_alt, nx, ny)
        az_x, az_y = get_line_points(m_az, c_az, nx, ny)

        # Plot the lines with clear labels
        ax.plot(alt_x, alt_y, 'r-', linewidth=2, label='Azimuth Line')
        ax.plot(az_x, az_y, 'b-', linewidth=2, label='Altitude Line')

        # Mark reference pixel
        ax.scatter(ref_x, ref_y, color='green', marker='x', s=200, label='Optical axis')

        # Labels and legend
        ax.set_title(f"Alt/Az Grid\nRot Angle: {rot_angle:.2f}", fontsize=16)
        ax.legend(loc='upper right', fontsize=12, frameon=True, facecolor='white', edgecolor='black')

        ax.text(0.05, 0.95, "N ↑  E ←", transform=ax.transAxes,
                fontsize=14, color='gold', fontweight='bold',
                ha='left', va='top', bbox=dict(facecolor='black', edgecolor='white', alpha=0.5))

        # Set axis limits to exact image dimensions
        ax.set_xlim(0, nx)
        ax.set_ylim(0, ny)

        # Save the output image with excellent resolution
        save_path = f"{output_dir}/{solved_fits.split('/')[-1]}.png"
        plt.savefig(save_path, dpi=400, bbox_inches='tight')
        plt.show()
        print(f"Saved image to {save_path}")

        # Subtract background for star detection
        mean, median, std = sigma_clipped_stats(image_data, sigma=5.0)
        image_data = image_data.astype(np.float64)
        image_data -= median  # Background subtraction

        # Detect stars using DAOStarFinder
        daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold * std)
        sources = daofind(image_data)
        print(f"Detected {len(sources)} sources")

        if sources is None or len(sources) < num_stars_to_analyze:
            print("Warning: Not enough stars detected in the image!")
            return

        # Select random stars
        np.random.seed(random_seed)  # For reproducibility
        selected_indices = np.random.choice(len(sources), size=num_stars_to_analyze, replace=False)
        selected_stars = sources[selected_indices]

        # Extract filename
        filename = solved_fits.split('/')[-1]

        # Create a figure with 4x3 grid for the subplots (changed from 2x2 to 4x3 for 12 stars)
        fig = plt.figure(figsize=(16, 16))  # Increased figure size to accommodate more subplots
        fig.suptitle(f"Stars Elongation Analysis - {filename} \nAlt-Az frame Rot Angle: {rot_angle:.2f}", fontsize=14)

        # Create subplots one by one with their own WCS
        for i, star in enumerate(selected_stars):
            # Extract star centroid
            x, y = star['xcentroid'], star['ycentroid']

            # Crop the image around the star and subtract background
            xmin = int(max(0, x - size))
            xmax = int(min(nx, x + size))
            ymin = int(max(0, y - size))
            ymax = int(min(ny, y + size))
            cropped_data = image_data[ymin:ymax, xmin:xmax]
            cropped_mean, cropped_median, _ = sigma_clipped_stats(cropped_data, sigma=5.0)
            cropped_data -= cropped_median  # Local background subtraction

            # Create a proper WCS object for the cropped region
            cropped_wcs = wcs.deepcopy()
            cropped_wcs.wcs.crpix[0] -= xmin
            cropped_wcs.wcs.crpix[1] -= ymin

            # Create a subplot with the proper WCS projection (4x3 grid)
            ax = fig.add_subplot(4, 3, i+1, projection=cropped_wcs)
            
            # Compute morphological properties
            cat = data_properties(cropped_data)

            # Extract parameters
            semimajor = cat.semimajor_sigma.value * r
            semiminor = cat.semiminor_sigma.value * r
            theta = cat.orientation.to(u.deg).value
            elongation = 1 - (semiminor / semimajor)

            # Display the cropped image
            zscale = ZScaleInterval()
            vmin, vmax = zscale.get_limits(cropped_data)
            ax.imshow(cropped_data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)

            # Correct ellipse position using offsets
            xypos = (cat.xcentroid, cat.ycentroid)  # Centroid in cropped frame
            aperture = EllipticalAperture(xypos, semimajor, semiminor, theta=cat.orientation.to(u.rad).value)
            aperture.plot(ax=ax, color='red', lw=2)

            # Compute and plot altitude and azimuth lines
            m_alt = np.tan(rot_angle.to(u.rad).value)
            m_az = -1 / m_alt
            # Calculate c based on the global x,y but then adjust for the cropped frame
            c_alt = y - m_alt * x  # Global intercept
            c_az = y - m_az * x    # Global intercept
            # Transform to cropped frame intercept
            c_alt_local = c_alt - ymin + m_alt * xmin
            c_az_local = c_az - ymin + m_az * xmin

            def get_line_points(m, c, width, height):
                points = []
                for x_edge in [0, width-1]:
                    y_edge = m * x_edge + c
                    if 0 <= y_edge < height:
                        points.append((x_edge, y_edge))
                for y_edge in [0, height-1]:
                    x_edge = (y_edge - c) / m if m != 0 else width/2
                    if 0 <= x_edge < width:
                        points.append((x_edge, y_edge))
                return np.array(points).T

            # Get cropped image dimensions
            crop_height, crop_width = cropped_data.shape
            alt_x, alt_y = get_line_points(m_alt, c_alt_local, crop_width, crop_height)
            az_x, az_y = get_line_points(m_az, c_az_local, crop_width, crop_height)
            
            # Plot the lines
            if len(alt_x) >= 2:
                ax.plot(alt_x, alt_y, 'r-', linewidth=2, label='Azimuth Line')
            if len(az_x) >= 2:
                ax.plot(az_x, az_y, 'b-', linewidth=2, label='Altitude Line')

            # Add North-East directions at **top-left**
            ax.text(0.05, 0.95, "N ↑  E ←", transform=ax.transAxes,
                    fontsize=12, color='yellow', fontweight='bold',
                    ha='left', va='top', bbox=dict(facecolor='black', edgecolor='white', alpha=0.5))

            # Add elongation and angle info at **bottom-left**
            ax.text(0.05, 0.05, f"Ellipse elongation: {elongation:.2f}\nEllipse angle: {theta:.2f}°",
                    transform=ax.transAxes, fontsize=8, color='white',
                    ha='left', va='bottom', bbox=dict(facecolor='black', edgecolor='white', alpha=0.5))

            # Add celestial grid
            ax.coords.grid(color='yellow', linestyle='dotted', linewidth=1)
            
            # Get RA/Dec of star center
            ra, dec = cropped_wcs.wcs_pix2world(cat.xcentroid, cat.ycentroid, 0)

            ax.text(0.95, 0.05, f"Star's RA: {ra:.4f}\n Star's Dec: {dec:.4f}", transform=ax.transAxes,
                fontsize=8, color='white', ha='right', va='bottom',
                bbox=dict(facecolor='black', edgecolor='white', alpha=0.5))
            
            # Set axis limits to the cropped region
            ax.set_xlim(0, crop_width*0.98)
            ax.set_ylim(0, crop_height*0.98)


            # Add legend for azimuth and altitude lines
            ax.legend(loc="upper right", fontsize=8, facecolor='white', edgecolor='black', framealpha=1)

        # Adjust layout and save the figure
        plt.tight_layout()
        save_path = f"{output_dir}/{filename}_stars.png"
        plt.savefig(save_path, dpi=400, bbox_inches='tight')
        plt.show()
        print(f"Saved image to {save_path}")

# Example usage
plot_with_grids_and_stars(input_fits, output_dir, num_stars_to_analyze, fwhm, threshold, random_seed)