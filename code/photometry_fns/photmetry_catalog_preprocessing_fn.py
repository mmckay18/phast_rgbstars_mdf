# pre_processing scripts for processing CMD isochrone tables

import pandas as pd
import os
import re
from tqdm import tqdm
import glob
import numpy as np
from scipy.interpolate import LinearNDInterpolator
from astropy.io import fits
from astropy.table import Table
import vaex


def make_optical_photmetry_from_hst_dataproduct(
    photmetry_catalog_filepath, output_csv_filepath
):
    """
    Takes the HST photmetry catalog and filters it down to just the optical bands for RGB metallicity istrubution
    analysis. The output is saved as a csv.

    """
    if photmetry_catalog_filepath.endswith(".fits"):

        phast_catalog = fits.open(photmetry_catalog_filepath)

        # phast_catalog.info()
        # phast_catalog[1].header['*']

        phast_table = phast_catalog[1].data
        # phast_table.shape

        print("Adding F475W and F814W columns to table")
        phast_table = Table()
        phast_table["ra"] = phast_catalog[1].data["ra"]
        phast_table["dec"] = phast_catalog[1].data["dec"]

        # F475W
        phast_table["f475w_vega"] = phast_catalog[1].data["f475w_vega"]
        phast_table["f475w_snr"] = phast_catalog[1].data["f475w_snr"]
        phast_table["f475w_crowd"] = phast_catalog[1].data["f475w_crowd"]
        phast_table["f475w_sharp"] = phast_catalog[1].data["f475w_sharp"]
        phast_table["f475w_flag"] = phast_catalog[1].data["f475w_flag"]

        # F814W
        phast_table["f814w_vega"] = phast_catalog[1].data["f814w_vega"]
        phast_table["f814w_snr"] = phast_catalog[1].data["f814w_snr"]
        phast_table["f814w_crowd"] = phast_catalog[1].data["f814w_crowd"]
        phast_table["f814w_sharp"] = phast_catalog[1].data["f814w_sharp"]
        phast_table["f814w_flag"] = phast_catalog[1].data["f814w_flag"]

        # F275W
        phast_table["f275w_vega"] = phast_catalog[1].data["f275w_vega"]
        phast_table["f275w_snr"] = phast_catalog[1].data["f275w_snr"]
        phast_table["f275w_crowd"] = phast_catalog[1].data["f275w_crowd"]
        phast_table["f275w_sharp"] = phast_catalog[1].data["f275w_sharp"]
        phast_table["f275w_flag"] = phast_catalog[1].data["f275w_flag"]

        print("Save catalog to CSV")
        df = phast_table.to_pandas()
        df.to_csv(f"{output_csv_filepath}", index=True)

    elif photmetry_catalog_filepath.endswith(".hdf5"):
        ds = vaex.open(photmetry_catalog_filepath)
        (
            ra,
            dec,
            f475w_vega,
            f475w_snr,
            f475w_crowd,
            f475w_sharp,
            f475w_flag,
            f814w_vega,
            f814w_snr,
            f814w_crowd,
            f814w_sharp,
            f814w_flag,
            f275w_vega,
            f275w_snr,
            f275w_crowd,
            f275w_sharp,
            f275w_flag,
        ) = ds.evaluate(
            [
                "RA",
                "DEC",
                "F475W_VEGA",
                "F475W_SNR",
                "F475W_CROWD",
                "F475W_SHARP",
                "F475W_GST_FLAG",
                "F814W_VEGA",
                "F814W_SNR",
                "F814W_CROWD",
                "F814W_SHARP",
                "F814W_GST_FLAG",
                "F275W_VEGA",
                "F275W_SNR",
                "F275W_CROWD",
                "F275W_SHARP",
                "F275W_GST_FLAG",
            ],
            # selection='(F475W_GST_FLAG == 0.0) & (F814W_GST_FLAG == 0.0) & (F475W_CROWD+F814W_CROWD <= 1.0) & (F475W_SNR >= 4.0) & (F814W_SNR >= 4.0)')
            selection="((F475W_SHARP+F814W_SHARP)**2 <= 0.075) & (F475W_CROWD+F814W_CROWD <= 1.0) & (F475W_SNR >= 4.0) & (F814W_SNR >= 4.0)",
        )

        ds.close()  # if you don't need anything else from the file

        # Save ds as a CSV file
        print("Writing to dataframe...")
        df = pd.DataFrame(
            {
                "ra": ra,
                "dec": dec,
                "f475w_vega": f475w_vega,
                "f475w_snr": f475w_snr,
                "f475w_crowd": f475w_crowd,
                "f475w_sharp": f475w_sharp,
                "f475w_flag": f475w_flag,
                "f814w_vega": f814w_vega,
                "f814w_snr": f814w_snr,
                "f814w_crowd": f814w_crowd,
                "f814w_sharp": f814w_sharp,
                "f814w_flag": f814w_flag,
                "f275w_vega": f275w_vega,
                "f275w_snr": f275w_snr,
                "f275w_crowd": f275w_crowd,
                "f275w_sharp": f275w_sharp,
                "f275w_flag": f275w_flag,
            }
        )
        # df = df.loc[((df['f475w_sharp']+df['f814w_sharp'])**2 <= 0.075)]
        print("Saving as CSV file...")
        df.to_csv(
            "/astro/users/mmckay18/analysis_routine_PHAST/DATA/phat_f475W_f814W_table.csv",
            index=False,
        )

        print(f"Total number of stars: {len(df)}")


def reduce_optical_photmetry_table(photometry_fits_filepath, output_filepath):
    """
    Reads in the photometry fits table produced by the HST pipeline and filters it using the criteria from
    Gregersen et al 2015 for GST. The remaining stars are foreground extinction corrected and saved as a csv.

    Parameters
    ----------
    photometry_fits_filepath : str
        Path to the photometry fits table produced by the HST pipeline
    output_filepath : str
        Path to the directory where the output csv will be saved
        (i.e /Users/mmckay/phd_projects/analysis_routine/DATA/reduced_phast_catalog.csv)

    Returns
    -------
    catalog_df : pandas dataframe
        Dataframe of the reduced catalog

    """
    df = pd.read_csv(photometry_fits_filepath)

    # Gregersen et al 2015 criteria GST(good stars)
    catalog_df = df.loc[
        (df["f475w_flag"] == 0.0)
        & (df["f814w_flag"] == 0.0)
        & (df["f475w_snr"] >= 4.0)
        & (df["f814w_snr"] >= 4.0)
        & (df["f475w_crowd"] + df["f814w_crowd"] <= 1.0)
        & ((df["f475w_sharp"] + df["f814w_sharp"]) ** 2 <= 0.075)
        #  & (df['f814w_vega'] <= 23.0) # Apply after foreground reddening correction
    ]

    # Apply foreground reddening correction
    # A_lam/A_V: F814W=0.59696, F475W=1.21182 from CMD 3.7 output A_V = 0.17
    catalog_df["f814w_vega_ecorr"] = catalog_df["f814w_vega"] + 0.596 * 0.17
    catalog_df["f475w_vega_ecorr"] = catalog_df["f475w_vega"] + 1.212 * 0.17
    catalog_df["f475w-f814w_ecorr"] = (
        catalog_df["f475w_vega_ecorr"] - catalog_df["f814w_vega_ecorr"]
    )
    catalog_df = catalog_df.loc[
        (catalog_df["f814w_vega_ecorr"] <= 23.0)
    ]  # Apply after foreground reddening correction

    # Save the table
    catalog_df.to_csv(output_filepath)

    return catalog_df


def phast_reduce_optical_photmetry_table(photometry_fits_filepath, output_filepath):
    """
    Reads in the photometry fits table produced by the HST pipeline and filters it using the criteria from
    Gregersen et al 2015 for GST. The remaining stars are foreground extinction corrected and saved as a csv.

    Parameters
    ----------
    photometry_fits_filepath : str
        Path to the photometry fits table produced by the HST pipeline
    output_filepath : str
        Path to the directory where the output csv will be saved
        (i.e /Users/mmckay/phd_projects/analysis_routine/DATA/reduced_phast_catalog.csv)

    Returns
    -------
    catalog_df : pandas dataframe
        Dataframe of the reduced catalog

    """
    df = pd.read_csv(photometry_fits_filepath)

    # Gregersen et al 2015 criteria GST(good stars)
    catalog_df = df.loc[
        (df["f475w_flag"] < 8.0) &
        (df["f814w_flag"] < 8.0)
        & (df["f475w_snr"] >= 4.0)
        & (df["f814w_snr"] >= 4.0)
        & (df["f475w_crowd"] + df["f814w_crowd"] <= 1.0)
        & ((df["f475w_sharp"] + df["f814w_sharp"]) ** 2 <= 0.075)
        #  & (df['f814w_vega'] <= 23.0) # Apply after foreground reddening correction
    ]

    # Apply foreground reddening correction
    # A_lam/A_V: F814W=0.59696, F475W=1.21182 from CMD 3.7 output A_V = 0.17
    catalog_df["f814w_vega_ecorr"] = catalog_df["f814w_vega"] + 0.596 * 0.17
    catalog_df["f475w_vega_ecorr"] = catalog_df["f475w_vega"] + 1.212 * 0.17
    catalog_df["f475w-f814w_ecorr"] = (
        catalog_df["f475w_vega_ecorr"] - catalog_df["f814w_vega_ecorr"]
    )
    catalog_df = catalog_df.loc[
        (catalog_df["f814w_vega_ecorr"] <= 23.0)
    ]  # Apply after foreground reddening correction

    # Save the table
    catalog_df.to_csv(output_filepath)

    return catalog_df


def catalog_linear_interpolation(
    catalog_csv_filepath, isochrone_csv_filepath, output_filepath, savefile=True
):
    """
    Interpolates the isochrone table to the catalog table using the F814W_appmag and F475W_appmag-F814W_appmag
    columns from the isochrone table and the f814w_vega and f475w_vega-f814w_vega columns from the catalog table.
    The interpolated values are then added as new columns to the catalog table and saved as a new csv.

    Parameters
    ----------
    catalog_csv_filepath : str
        Path to the catalog csv produced by make_optical_photmetry_table
    isochrone_csv_filepath : str
        Path to the isochrone csv produced by read_sort_isochrone_table
    output_filepath : str
        Path to the directory where the output csv will be saved
        (i.e /Users/mmckay/phd_projects/analysis_routine/DATA/interpolated_phast_MH_catalog.csv)
    savefile : bool, optional
        If True, saves the catalog with interpolated values as a csv, by default True


    Returns
    -------
    catalog_df : pandas dataframe
        Dataframe of all catalog stars with interpolated isochrone values
    """
    # Read in the catalog and isochrone tables
    catalog_df = pd.read_csv(catalog_csv_filepath)
    RGB_merged_df = pd.read_csv(isochrone_csv_filepath)

    # Create a 2D array of F814W and F475W-F814W magnitudes from RGB_merged_df
    cmd_arr = RGB_merged_df[["F814W_appmag", "F475W_appmag-F814W_appmag"]].values

    # Mini,int_IMF,Mass,logL,logTe,logg,
    # Create a LinearNDInterpolator object with 'MH' as the values
    print("Running LinearNDInterpolator")
    interp_MH = LinearNDInterpolator(cmd_arr, RGB_merged_df["MH"], fill_value=-99.0)
    interp_Mini = LinearNDInterpolator(cmd_arr, RGB_merged_df["Mini"], fill_value=-99.0)
    interp_int_IMF = LinearNDInterpolator(
        cmd_arr, RGB_merged_df["int_IMF"], fill_value=-99.0
    )
    interp_Mass = LinearNDInterpolator(cmd_arr, RGB_merged_df["Mass"], fill_value=-99.0)
    interp_logL = LinearNDInterpolator(cmd_arr, RGB_merged_df["logL"], fill_value=-99.0)
    interp_logTe = LinearNDInterpolator(
        cmd_arr, RGB_merged_df["logTe"], fill_value=-99.0
    )
    interp_logg = LinearNDInterpolator(cmd_arr, RGB_merged_df["logg"], fill_value=-99.0)

    # Create a 2D array of f814w_vega and f475w-f814w magnitudes from catalog_df
    # catalog_df['f475w_vega-f814w_vega'] = catalog_df['f475w_vega'] - catalog_df['f814w_vega']

    # PHAST catalog color magnitude diagram values
    filtered_arr = catalog_df[["f814w_vega_ecorr", "f475w-f814w_ecorr"]].values

    # Interpolate the 'MH' values for each row in catalog_df
    interpolated_MH = interp_MH(filtered_arr)
    interpolated_Mini = interp_Mini(filtered_arr)
    interpolated_int_IMF = interp_int_IMF(filtered_arr)
    interpolated_Mass = interp_Mass(filtered_arr)
    interpolated_logL = interp_logL(filtered_arr)
    interpolated_logTe = interp_logTe(filtered_arr)
    interpolated_logg = interp_logg(filtered_arr)

    # Add the interpolated 'MH' values as a new column in catalog_df
    catalog_df["interpolated_MH"] = interpolated_MH
    catalog_df["interpolated_Mini"] = interpolated_Mini
    catalog_df["interpolated_int_IMF"] = interpolated_int_IMF
    catalog_df["interpolated_Mass"] = interpolated_Mass
    catalog_df["interpolated_logL"] = interpolated_logL
    catalog_df["interpolated_logTe"] = interpolated_logTe
    catalog_df["interpolated_logg"] = interpolated_logg

    print(
        f'Dropping fill values if any {len(catalog_df[catalog_df["interpolated_MH"] == -99.0])}'
    )
    catalog_df = catalog_df[catalog_df["interpolated_MH"] != -99.0].dropna()

    # Drop  M/H values > 0.6 ~ CHECK IF THIS THE RIGHT THING TO DO
    # catalog_df = catalog_df[catalog_df['interpolated_MH'] <= 0.6].dropna()

    # Save the merged dataframe as CSV
    if savefile is True:
        print(
            f"Saving the photometry catalog with interpolated values as {output_filepath}"
        )
        catalog_df.to_csv(f"{output_filepath}", index=False)
    else:
        pass

    return catalog_df
