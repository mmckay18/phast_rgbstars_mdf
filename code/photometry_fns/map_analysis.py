from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os


os.chdir('/Users/mmckay/phd_projects/analysis_routine/code/photometry_fns')
from plotting_fns import (create_spatial_bins_and_median_optimized)


def align_images(input_file_high_res, input_file_low_res):
    # Open the high-resolution FITS file
    hdul_high_res = fits.open(input_file_high_res)
    data_high_res = hdul_high_res[0].data.T
    wcs_high_res = WCS(hdul_high_res[0].header)

    # Open the low-resolution FITS file
    hdul_low_res = fits.open(input_file_low_res)
    data_low_res = hdul_low_res[0].data
    wcs_low_res = WCS(hdul_low_res[0].header)

    # Reproject the low-resolution image to match the high-resolution WCS
    data_reprojected, _ = reproject_interp((data_low_res, wcs_low_res), wcs_high_res, shape_out=data_high_res.shape)

    # Close the FITS files
    hdul_high_res.close()
    hdul_low_res.close()

    return data_high_res, data_reprojected, wcs_high_res, data_low_res


def set_phast_wcs(input_FITS_filepath, overwrite=False, show_plot=True):
    phast_hdu = fits.open('/Users/mmckay/phd_projects/analysis_routine/DATA/phast_cmd_mdf_spatial_gregersen_box.fits')

    # astropy create for binned maps
    phast_hdu[0].header['*']
    # reverse x column to match dust image
    # phast_hdu[0].data = np.fliplr(phast_hdu[0].data)
    w = WCS(phast_hdu[0].header)
    w.wcs.crpix = [80, 92] # center pixel
    # w.wcs.crpix = [0, 0] # center pixel
    w.wcs.crval = [10.79, 41.27] # RA and dec values in hours and degreesixel scale in degrees/pixel
    # w.wcs.crval = [10.79, 41.13] # RA and dec values in hours and degreesixel scale in degrees/pixel
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"] # coordinate system
    # w.wcs.cdelt = [-0.005, 0.005] # pixel scale in degrees/pixel
    w.wcs.cdelt = [0.01, 0.01] # pixel scale in degrees/pixel
    print(w)

    phast_hdu[0].header.update(w.to_header())
    phast_hdu.writeto('input_FITS_filepath', overwrite=overwrite)

    if show_plot:
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection=w)

        # Plot the image data
        ax.imshow(phast_hdu[0].data, cmap='magma', origin='lower')
        # Add contour lines
        # contour_levels = np.nanpercentile(phast_hdu[0].data, [-1])
        # plt.contour(phast_hdu[0].data, levels=contour_levels, colors='black')

        # ax.imshow(dust_lum_hdu[0].data, cmap='magma', origin='lower')
        # Set the axis labels
        ax.set_xlabel('RA')
        ax.set_ylabel('Dec')

        # Show the plot
        plt.show()
        phast_hdu.close()

    else:
        phast_hdu.close()


def set_m31_wcs(input_FITS_filepath, overwrite=False, show_plot=True):
    m31_hdu = fits.open(input_FITS_filepath, mode='update')

    #! Remove PC1_1 and PC2_2 from the header
    if 'PC1_1' in m31_hdu[0].header:
        m31_hdu[0].header.remove('PC1_1')
    if 'PC1_2' in m31_hdu[0].header:
        m31_hdu[0].header.remove('PC1_2')
    if 'PC2_1' in m31_hdu[0].header:
        m31_hdu[0].header.remove('PC2_1')
    if 'PC2_2' in m31_hdu[0].header:
        m31_hdu[0].header.remove('PC2_2')

    if 'CD1_1' in m31_hdu[0].header:
        m31_hdu[0].header.remove('CD1_1')
    if 'CD1_2' in m31_hdu[0].header:
        m31_hdu[0].header.remove('CD1_2')
    if 'CD2_1' in m31_hdu[0].header:
        m31_hdu[0].header.remove('CD2_1')
    if 'CD2_2' in m31_hdu[0].header:
        m31_hdu[0].header.remove('CD2_2')


    # astropy create for binned maps
    m31_hdu[0].header = m31_hdu[0].header['*'] # reset header to default

    # m31_hdu[0].header['*']
    # reverse x column to match dust image
    # m31_hdu[0].data = np.fliplr(m31_hdu[0].data)
    w = WCS(m31_hdu[0].header)
    # w.wcs.crpix = [93, 70] # center pixel - align with PHAT dust structure (old)
    w.wcs.crpix = [84, 81] # center pixel - align with PHAT dust structure
    w.wcs.crval = [10.79, 41.13] # RA and dec values in hours and degreesixel scale in degrees/pixel
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"] # coordinate system
    w.wcs.cdelt = [np.rad2deg(np.radians(0.01)*np.cos(np.radians(w.wcs.crval[1]))), 0.01] # pixel scale in degrees/pixel #* Works with the right scaling...
    print(w)

    m31_hdu[0].header.update(w.to_header())
    m31_hdu.writeto(input_FITS_filepath, overwrite=overwrite)

    if show_plot:
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection=w)

        # Plot the image data
        ax.imshow(m31_hdu[1].data, cmap='magma', origin='lower')
        # Add contour lines
        # contour_levels = np.nanpercentile(m31_hdu[0].data, [-1])
        # plt.contour(m31_hdu[0].data, levels=contour_levels, colors='black')

        # ax.imshow(dust_lum_hdu[0].data, cmap='magma', origin='lower')
        # Set the axis labels
        ax.set_xlabel('RA')
        ax.set_ylabel('Dec')

        # Show the plot
        plt.show()
        m31_hdu.close()

    else:
        m31_hdu.close()


def make_m31_maps_FITS(csv_filepath, output_fitsfile='/Users/mmckay/phd_projects/analysis_routine/DATA/all_catalog_maps_fitsfiles/test_m31_catalog.fits', overwrite=False):
    """
    Create FITS maps from a CSV file containing photometry catalog, isochrone estimates.

    Parameters:
    - csv_filepath (str): The file path to the CSV file containing the photometry catalog and metallicity estimates.
    - output_fitsfile (str): The file path to save the output FITS file. Default is '/Users/mmckay/phd_projects/analysis_routine/DATA/all_catalog_maps_fitsfiles/test_m31_catalog.fits'.
    - overwrite (bool): Whether to overwrite the output FITS file if it already exists. Default is False.

    Returns:
    None
    

    Example:
    from map_analysis import *
    csv_filepath = '/Users/mmckay/phd_projects/analysis_routine/DATA/interpolated_m31_MH_catalog.csv'
    make_m31_maps_FITS(csv_filepath, output_fitsfile='/Users/mmckay/phd_projects/analysis_routine/DATA/all_catalog_maps_fitsfiles/test_m31_catalog.fits', overwrite=True)
    """
    
    # Read in the dataframe
    print('Reading catalog...')
    df = pd.read_csv(csv_filepath)
    
    # Create a new FITS file
    hdu = fits.PrimaryHDU()
    hdul = fits.HDUList([hdu])
    # Loop over each column in the dataframe
    print('Binning Maps...')
    for col_name in df.columns:
        if col_name not in ['ra', 'dec', 'Unnamed: 0']:
            print(col_name.upper())
            ra = df['ra']
            dec = df['dec']
            col_data = df[col_name]
            try:
                print('Binning...')
                min_ra, max_ra, min_dec, max_dec, median_map, density_map, bin_counts_nonzero = create_spatial_bins_and_median_optimized(ra, dec, col_data, bin_size_deg=0.01)
                # Create a new image extension
                hdu_new = fits.ImageHDU(median_map.T)
                print(median_map.shape)
                # Set the header information
                hdu_new.header['EXTNAME'] = col_name.upper()
                # hdu_new.header['COMMENT'] = 'Median map for column: {}'.format(col)
                # Append the new extension to the HDU list
                print('Appending Binned Map to HDU object')
                hdul.append(hdu_new)
            except Exception as e:
                print(f'Error with column: {col_name}. Exception: {e}')
                continue
        else:
            pass


    #Save fits file
    print('Saving fits file...')
    hdul.writeto(output_fitsfile, overwrite=overwrite)
    print(f'Savings {output_fitsfile}')
    hdul.close()


def calculate_elliptical_radius(mh_map, center_x, center_y, semi_major_axis, semi_minor_axis, angle):
    """
    Calculates the elliptical radius for each point in a map given the center coordinates,
    semi-major axis, semi-minor axis, and rotation angle of the ellipse.
    
    Parameters:
    - mh_map (numpy.ndarray): 2D numpy array representing the map on which to calculate the elliptical radius.
    - center_x (float): The x-coordinate of the center of the ellipse.
    - center_y (float): The y-coordinate of the center of the ellipse.
    - semi_major_axis (float): Length of the semi-major axis of the ellipse.
    - semi_minor_axis (float): Length of the semi-minor axis of the ellipse.
    - angle (float): Rotation angle of the ellipse in radians.
    
    Returns:
    - numpy.ndarray: 2D numpy array of the same shape as mh_map, where each element
      represents the elliptical radius of that point from the center of the ellipse.
    """
    import numpy as np
    
    # Dimensions of the mh_map
    height, width = mh_map.shape

    # Create a grid of x, y coordinates
    x = np.arange(width)
    y = np.arange(height)
    x_grid, y_grid = np.meshgrid(x, y)

    # Translate grid to be centered at (center_x, center_y)
    x_grid_centered = x_grid - center_x
    y_grid_centered = y_grid - center_y

    # Rotate the grid by the angle to align with the major and minor axes
    cos_angle = np.cos(angle)
    sin_angle = np.sin(angle)
    x_rotated = cos_angle * x_grid_centered + sin_angle * y_grid_centered
    y_rotated = -sin_angle * x_grid_centered + cos_angle * y_grid_centered

    # Calculate the elliptical radius for each point
    elliptical_radius_map = np.sqrt((x_rotated / semi_major_axis)**2 + (y_rotated / semi_minor_axis)**2)
    
    return elliptical_radius_map


def create_elliptical_mask(center_x, center_y, semi_major_axis, semi_minor_axis, angle, array_shape, inner_radius, outer_radius):
    """
    Creates a boolean mask for an elliptical annulus based on specified parameters.
    
    Parameters:
    - center_x (float): The x-coordinate of the center of the ellipse.
    - center_y (float): The y-coordinate of the center of the ellipse.
    - semi_major_axis (float): Length of the semi-major axis of the ellipse.
    - semi_minor_axis (float): Length of the semi-minor axis of the ellipse.
    - angle (float): Rotation angle of the ellipse in radians, measured from the x-axis.
    - array_shape (tuple of int): The shape of the array (height, width) for which the mask is created.
    - inner_radius (float): The inner radius of the elliptical annulus.
    - outer_radius (float): The outer radius of the elliptical annulus.
    
    Returns:
    - numpy.ndarray: A 2D boolean array where True represents points within the elliptical annulus.
    """
    import numpy as np
    
    # Generate a grid of coordinates centered at (center_x, center_y)
    y, x = np.ogrid[-center_y:array_shape[0]-center_y, -center_x:array_shape[1]-center_x]
    
    # Calculate the ellipse equation components
    ellipse = ((x * np.cos(angle) + y * np.sin(angle))**2 / semi_major_axis**2 +
               (y * np.cos(angle) - x * np.sin(angle))**2 / semi_minor_axis**2)
    
    # Create the mask based on the inner and outer radius conditions
    mask = (ellipse >= inner_radius**2) & (ellipse <= outer_radius**2)
    
    return mask