from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# FILEPATH: /Users/mmckay/phd_projects/analysis_routine/code/phat_phast_analysis.ipynb
def create_spatial_bins_and_median_optimized(ra, dec, values, bin_size_deg):
    # Convert RA and Dec to SkyCoord object
    coords = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)

    # Define the grid boundaries based on the region of interest
    min_ra, max_ra = min(coords.ra.deg), max(coords.ra.deg)
    min_dec, max_dec = min(coords.dec.deg), max(coords.dec.deg)

    # Calculate the number of bins in each dimension
    num_bins_ra = int((max_ra - min_ra) / bin_size_deg) + 1
    num_bins_dec = int((max_dec - min_dec) / bin_size_deg) + 1

    # Create bins using numpy.digitize
    ra_bins = np.digitize(coords.ra.deg, np.linspace(min_ra, max_ra, num_bins_ra))
    dec_bins = np.digitize(coords.dec.deg, np.linspace(min_dec, max_dec, num_bins_dec))

    # Calculate median values using numpy.bincount
    flat_indices = num_bins_dec * (ra_bins - 1) + dec_bins - 1
    bin_counts = np.bincount(flat_indices, minlength=num_bins_ra * num_bins_dec)
    bin_values = np.bincount(flat_indices, weights=values, minlength=num_bins_ra * num_bins_dec)

    # Avoid division by zero by setting zero counts to NaN
    bin_counts_nonzero = bin_counts.astype(float)
    bin_counts_nonzero[bin_counts_nonzero == 0] = np.nan

    # Calculate median values
    median_values = bin_values / bin_counts_nonzero
    median_values = median_values.reshape((num_bins_ra, num_bins_dec))

    return min_ra, max_ra, min_dec, max_dec, median_values

# Load CSV files
phat_df = pd.read_csv('/Users/mmckay/phd_projects/analysis_routine/DATA/interpolated_phat_MH_catalog.csv')
phast_df = pd.read_csv('/Users/mmckay/phd_projects/analysis_routine/DATA/interpolated_phast_MH_catalog.csv')


# PHAST
print("Defining example usage variables")
ra_data = phast_df['ra']  # Replace with your RA data
dec_data = phast_df['dec']  # Replace with your Dec data
values_data = phast_df['interpolated_MH']  # Replace with your values data corresponding to each object
bin_size_deg = 0.01  # Choose your bin size in degrees


print("Calculating spatial bins and median values")
min_ra, max_ra, min_dec, max_dec, median_values = create_spatial_bins_and_median_optimized(ra_data, dec_data, values_data, bin_size_deg)

print("Plotting the bins and color coding by median values")
# Plot the bins and color code by median values
plt.imshow(median_values.T, origin='lower', extent=[min_ra, max_ra, min_dec, max_dec], cmap='magma', label=f'PHAST {len(median_values)}')
plt.colorbar(label='Median Values')
plt.xlabel('RA')
plt.ylabel('Dec')
plt.title(f'Spatial Bin = {bin_size_deg}')
plt.legend()
plt.savefig('/Users/mmckay/phd_projects/analysis_routine/FIGURES/phast_MH_spatial_bins.jpg', format='jpeg')
plt.show()

# PHAT
ra_data = phat_df['ra']  # Replace with your RA data
dec_data = phat_df['dec']  # Replace with your Dec data
values_data = phat_df['interpolated_MH']  # Replace with your values data corresponding to each object
bin_size_deg = 0.01  # Choose your bin size in degrees

print("Defining example usage variables")
print("Calculating spatial bins and median values")
min_ra, max_ra, min_dec, max_dec, median_values = create_spatial_bins_and_median_optimized(ra_data, dec_data, values_data, bin_size_deg)


print("Plotting the bins and color coding by median values")# Plot the bins and color code by median values
plt.imshow(median_values.T, origin='lower', extent=[min_ra, max_ra, min_dec, max_dec], cmap='magma', label=f'PHAT {len(median_values)}')
plt.colorbar(label='Median Values')
plt.xlabel('RA')
plt.ylabel('Dec')
plt.title(f'Spatial Bin = {bin_size_deg}')
plt.legend()
plt.savefig('/Users/mmckay/phd_projects/analysis_routine/FIGURES/phat_MH_spatial_bins.jpg', format='jpeg')
plt.show()

