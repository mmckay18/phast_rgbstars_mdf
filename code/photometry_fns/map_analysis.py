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
    m31_hdu = fits.open('/Users/mmckay/phd_projects/analysis_routine/DATA/m31_cmd_mdf_spatial_gregersen_box.fits')

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
    w.wcs.crpix = [93, 70] # center pixel - align with PHAT dust structure
    w.wcs.crval = [10.79, 41.13] # RA and dec values in hours and degreesixel scale in degrees/pixel
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"] # coordinate system
    w.wcs.cdelt = [np.rad2deg(np.radians(0.01)*np.cos(np.radians(w.wcs.crval[1]))), 0.01] # pixel scale in degrees/pixel #* Works with the right scaling...
    print(w)

    m31_hdu[0].header.update(w.to_header())
    m31_hdu.writeto('input_FITS_filepath', overwrite=overwrite)

    if show_plot:
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection=w)

        # Plot the image data
        ax.imshow(m31_hdu[0].data, cmap='magma', origin='lower')
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
    """
    
    # Read in the dataframe
    print('Reading catalog...')
    df = pd.read_csv(csv_filepath)
    
    # Create a new FITS file
    hdu = fits.PrimaryHDU()
    hdul = fits.HDUList([hdu])
    # Loop over each column in the dataframe
    print('Binning Maps...')
    for col_name in df.columns[3:4]:
        if col_name not in ['ra', 'dec']:
            print(col_name.upper())
            ra = df['ra']
            dec = df['dec']
            col_data = df[col_name]
            try:
                print('Binning...')
                min_ra, max_ra, min_dec, max_dec, median_map, density_map, bin_counts_nonzero = create_spatial_bins_and_median_optimized(ra, dec, col_data, bin_size_deg=0.01)
                # Create a new image extension
                print(median_map.shape)
                hdu_new = fits.ImageHDU(median_map)
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
