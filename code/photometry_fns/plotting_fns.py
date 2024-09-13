import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import numpy as np
from scipy.special import expit, logit
from scipy.interpolate import interp1d
from astropy.coordinates import SkyCoord
import astropy.units as u
import os
from astropy.io import fits
import matplotlib.colors as colors


def plot_CMD(
    csv_filepath,
    f814w_col,
    f475w_f814_col,
    z_col,
    savefile=True,
    apply_phat_rgb_MH_selection=False,
    save_rgb_selection=False,
    save_filepath="/Users/mmckay/phd_projects/analysis_routine/DATA/rgb_slection.csv",
    savefig_filepath="/Users/mmckay/phd_projects/analysis_routine/DATA/cmd_plot.png",
):
    """
    Plots the color magnitude values for the optical bands

    """
    df = pd.read_csv(csv_filepath)
    # Show the plot
    # plt.show()

    # TODO: write function make artifical star test
    if apply_phat_rgb_MH_selection is True:
        # Gregersen et al 2015 criteria for RGB stars were cut using the artificial star test
        # CMD coordinates of boundaries from 2015 paper
        fig, ((ax0, ax1)) = plt.subplots(
            figsize=(18, 10), nrows=1, ncols=2, sharey=True
        )
        scatter = ax0.scatter(
            df[f475w_f814_col], df[f814w_col], c=df[z_col], cmap="YlOrBr", s=1
        )
        ax0.set_xlabel(f475w_f814_col)
        ax0.set_ylabel(f814w_col)
        ax0.set_ylim(23.5, 20)
        cbar = plt.colorbar(scatter)
        cbar.set_label(z_col)

        # Set outer boundary
        x = [1.5, 2.8, 4, 5, 6, 7]
        y = [23, 23, 21.75, 21.4, 21.25, 21.1]
        # Create an interpolation function
        interp_func = interp1d(x, y, bounds_error=False, fill_value="extrapolate")
        # Interpolate x and y values
        x_interp = np.linspace(min(x), max(x), 1000)
        y_interp = interp_func(x_interp)

        # Set inner boundary
        x1 = [1.5, 1.6, 1.8]
        y1 = [23, 21.7, 20.3]
        interp_func1 = interp1d(x1, y1, bounds_error=False, fill_value="extrapolate")
        x1_interp = np.linspace(min(x1), max(x1), 1000)
        y1_interp = interp_func1(x1_interp)

        x2 = [7, 6, 5, 4, 2.8, 1.8]
        y2 = [21.1, 20.7, 20.25, 20.15, 20.1, 20.3]
        interp_func2 = interp1d(x2, y2, bounds_error=False, fill_value="extrapolate")
        x2_interp = np.linspace(min(x2), max(x2), 1000)
        y2_interp = interp_func2(x2_interp)

        # Exclude values below the interpolated line
        excluded_df = df[
            (df[f814w_col] > interp_func(df[f475w_f814_col]))
            | (df[f814w_col] <= interp_func1(df[f475w_f814_col]))
            | (df[f814w_col] <= interp_func2(df[f475w_f814_col]))
        ]
        ax0.scatter(excluded_df[f475w_f814_col], excluded_df[f814w_col], color="black")

        # Keep points within selection
        df = df[df[f814w_col] < interp_func(df[f475w_f814_col])]
        df = df[df[f814w_col] > interp_func1(df[f475w_f814_col])]
        df = df[df[f814w_col] > interp_func2(df[f475w_f814_col])]
        # ax0.scatter(excluded_df[f475w_f814_col], excluded_df[f814w_col], color='black')
        ax0.plot(x_interp, y_interp, color="black", linestyle="-", linewidth=1)
        ax0.plot(x1_interp, y1_interp, color="black", linestyle="-", linewidth=1)
        ax0.plot(x2_interp, y2_interp, color="black", linestyle="-", linewidth=1)

        # ax1
        scatter = ax1.scatter(
            df[f475w_f814_col], df[f814w_col], c=df[z_col], cmap="YlOrBr", s=1
        )
        ax1.plot(x_interp, y_interp, color="black", linestyle="-", linewidth=1)
        ax1.plot(x1_interp, y1_interp, color="black", linestyle="-", linewidth=1)
        ax1.plot(x2_interp, y2_interp, color="black", linestyle="-", linewidth=1)

        if save_rgb_selection is True:
            df.to_csv(save_filepath, index=False)
            print("Saved RGB selection to {}".format(save_filepath))
        if savefile is True:
            plt.savefig(savefig_filepath, dpi=300)
            print("Saved CMD plot to {}".format(savefig_filepath))

    else:
        # Create the scatter plot
        scatter = plt.scatter(
            df[f475w_f814_col], df[f814w_col], c=df[z_col], cmap="YlOrBr", s=1
        )

        # Set the plot title and labels
        # plt.title('Scatter Plot')
        plt.xlabel(f475w_f814_col)
        plt.ylabel(f814w_col)

        # Invert the y-axis
        plt.gca().invert_yaxis()

        # Add colorbar
        cbar = plt.colorbar(scatter)
        cbar.set_label(z_col)

        # Save the plot
        if savefile is True:
            plt.savefig(savefig_filepath, dpi=300)
    plt.show()
    plt.clf()


def create_spatial_bins_and_median_optimized(ra, dec, values, bin_size_deg):
    # Convert RA and Dec to SkyCoord object
    coords = SkyCoord(ra=ra, dec=dec, unit=u.degree)

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
    bin_values = np.bincount(
        flat_indices, weights=values, minlength=num_bins_ra * num_bins_dec
    )

    # Avoid division by zero by setting zero counts to NaN
    bin_counts_nonzero = bin_counts.astype(float)
    bin_counts_nonzero[bin_counts_nonzero == 0] = np.nan

    # Calculate median values
    median_values = bin_values / bin_counts_nonzero
    median_values = median_values.reshape((num_bins_ra, num_bins_dec))

    # calaculate the sum of the bin values - number density
    sum_values = bin_values / bin_counts_nonzero
    sum_values = bin_values.reshape((num_bins_ra, num_bins_dec))
    sum_values = np.where(sum_values == 0.0, np.nan, sum_values)

    # Reshape bin_counts to the shape of the grid
    bin_counts_nonzero = bin_counts_nonzero.reshape((num_bins_ra, num_bins_dec))

    return (
        min_ra,
        max_ra,
        min_dec,
        max_dec,
        median_values,
        sum_values,
        bin_counts_nonzero,
    )


def photometry_binned_spatial_map(
    catalog_csv_filepath, z_col, output_dir, catalog_name="phast"
):
    """ """
    df = pd.read_csv(catalog_csv_filepath)
    ra_data = df["ra"]  # Replace with your RA data
    dec_data = df["dec"]  # Replace with your Dec data
    values_data = df[
        z_col
    ]  # Replace with your values data corresponding to each object
    bin_size_deg = 0.01  # Choose your bin size in degrees

    # makes output dir if it doesnt exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    else:
        pass

    min_ra, max_ra, min_dec, max_dec, median_values, sum_values, bin_counts_nonzero = (
        create_spatial_bins_and_median_optimized(
            ra_data, dec_data, values_data, bin_size_deg
        )
    )

    # Plot the bins and color code by median values
    # print("Plotting the bins and color coding by median values")
    plt.figure(figsize=(18, 10))
    plt.imshow(
        median_values.T,
        origin="lower",
        extent=[min_ra, max_ra, min_dec, max_dec],
        cmap="magma",
    )
    plt.colorbar(label=f"Median {z_col}")
    plt.xlabel("RA")
    plt.ylabel("Dec")
    plt.title(f"Spatial Binning:{bin_size_deg}")
    plt.savefig(f"{output_dir}/{catalog_name}_{z_col}_median_binned_spatial_map.jpeg", dpi=300)
    plt.show()
    plt.clf()

    # Plot the bins and color code by sum values
    # print("Plotting the bins and color coding by sum values")
    plt.figure(figsize=(18, 10))
    plt.imshow(
        bin_counts_nonzero.T,
        origin="lower",
        extent=[min_ra, max_ra, min_dec, max_dec],
        cmap="ocean",
    )
    plt.colorbar(label=f"Sum {z_col}")
    plt.xlabel("RA")
    plt.ylabel("Dec")
    plt.title(f"Spatial Binning:{bin_size_deg}")
    plt.savefig(f"{output_dir}/{catalog_name}_{z_col}_sum_binned_spatial_map.jpeg", dpi=300)
    plt.show()
    plt.clf()


def plot_cmd_mdf_spatial(
    interpolated_phatter_rgb_agb,
    savefig_filepath,
    bin_size_deg=0.005,
    z_col="interpolated_MH",
):

    # Create a figure and subplots
    fig, axs = plt.subplots(1, 4, figsize=(25, 8))

    # Convert the object to numpy array before indexing
    f475w_f814w_ecorr = np.array(interpolated_phatter_rgb_agb["f475w-f814w_ecorr"])
    f814w_vega_ecorr = np.array(interpolated_phatter_rgb_agb["f814w_vega_ecorr"])
    z_col_arr = np.array(interpolated_phatter_rgb_agb[z_col])
    ra = np.array(interpolated_phatter_rgb_agb["ra"])
    dec = np.array(interpolated_phatter_rgb_agb["dec"])

    # Plot CMD
    # sns.scatterplot(x=f475w_f814w_ecorr, y=f814w_vega_ecorr, hue=z_col_arr, ax=axs[0], zorder=0, size=1, palette='magma')
    sns.kdeplot(
        x=f475w_f814w_ecorr,
        y=f814w_vega_ecorr,
        zorder=1,
        n_levels=10,
        thresh=0.01,
        fill=True,
        cbar=False,
        color="grey",
        ax=axs[0],
    )
    axs[0].set_xlabel("F475W-F814W", fontsize=14)
    axs[0].set_ylabel("F814W", fontsize=14)
    axs[0].set_title(f"Color-Magnitude Diagram N={len(f475w_f814w_ecorr)}", fontsize=14)
    axs[0].invert_yaxis()

    # Plot interpolated_MH KDE plot
    sns.histplot(data=z_col_arr, color="blue", kde=False, bins=10, ax=axs[1])
    axs[1].set_xlabel("M/H", fontsize=14)
    axs[1].set_ylabel("N", fontsize=14)
    axs[1].set_title("Metallicity Distribution Function (MDF)", fontsize=14)

    # Stellar Density binned map
    min_ra, max_ra, min_dec, max_dec, median_values, sum_values, density_map = (
        create_spatial_bins_and_median_optimized(
            ra, dec, z_col_arr, bin_size_deg=bin_size_deg
        )
    )
    im = axs[2].imshow(
        density_map.T,
        cmap="ocean",
        aspect="auto",
        extent=[min_ra, max_ra, max_dec, min_dec],
        norm=colors.LogNorm(),
    )
    axs[2].set_xlabel("RA", fontsize=14)
    axs[2].set_ylabel("DEC", fontsize=14)
    axs[2].set_title(f"Stellar Density ({bin_size_deg} deg^2)", fontsize=14)
    axs[2].invert_xaxis()
    axs[2].invert_yaxis()
    plt.colorbar(im, label=r"$\log(N \times 0.01\,\mathrm{deg}^2)$")

    # Median Metallicity Binned Spatial map
    im1 = axs[3].imshow(
        median_values.T,
        cmap="magma",
        aspect="auto",
        extent=[min_ra, max_ra, max_dec, min_dec],
    )
    axs[3].set_xlabel("RA", fontsize=14)
    axs[3].set_ylabel("DEC", fontsize=14)
    axs[3].set_title(f"Median M/H ({bin_size_deg} deg^2)", fontsize=14)
    axs[3].invert_xaxis()
    axs[3].invert_yaxis()
    plt.colorbar(im1, label=f"Median {z_col}")

    # Add a contour plot
    # contours = axs[3].contour(median_values.T, colors='black', origin='lower', levels=5)
    # plt.clabel(contours, inline=False, fontsize=8)
    # plt.colorbar(im).add_lines(contours)
    # Adjust spacing between subplots
    plt.tight_layout()
    plt.savefig(savefig_filepath, dpi=300)
    plt.show()

    # Save median M/H and stellar density maps to FITS file
    hdu_bin_counts = fits.PrimaryHDU(density_map.T)
    hdu_median = fits.ImageHDU(median_values.T)
    hdul = fits.HDUList([hdu_bin_counts, hdu_median])
    fitspath = savefig_filepath.split(".")[0] + ".fits"
    hdul.writeto(fitspath, overwrite=True)

    return median_values, density_map


def photometry_spatial_map():
    """ """


def mdf_histogram_plot():
    """ """
