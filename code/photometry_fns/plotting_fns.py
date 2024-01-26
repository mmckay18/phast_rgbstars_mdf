import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import numpy as np
from scipy.special import expit, logit
from scipy.interpolate import interp1d



def plot_CMD(csv_filepath, f814w_col, f475w_f814_col, z_col, savefile=True, apply_phat_rgb_MH_selection=False, save_rgb_selection=False, save_filepath='/Users/mmckay/phd_projects/analysis_routine/DATA/rgb_slection.csv', savefig_filepath='/Users/mmckay/phd_projects/analysis_routine/DATA/cmd_plot.png'):
    """
    Plots the color magnitude values for the optical bands
    
    """
    df = pd.read_csv(csv_filepath)
    # Show the plot
    # plt.show()
        
        # TODO: write function make artifical star test
    if apply_phat_rgb_MH_selection == True:
        # Gregersen et al 2015 criteria for RGB stars were cut using the artificial star test
        # CMD coordinates of boundaries from 2015 paper
        fig, ((ax0, ax1)) = plt.subplots(figsize=(18, 10), nrows=1, ncols=2, sharey=True)
        scatter = ax0.scatter(df[f475w_f814_col], df[f814w_col], c=df[z_col], cmap='YlOrBr', s=1)
        ax0.set_xlabel(f475w_f814_col)
        ax0.set_ylabel(f814w_col)
        ax0.set_ylim(23.5, 20)
        cbar = plt.colorbar(scatter)
        cbar.set_label(z_col)

        # Set outer boundary
        x = [1.5, 2.8, 4, 5, 6, 7]
        y = [23, 23, 21.75, 21.4, 21.25, 21.1]
        # Create an interpolation function
        interp_func = interp1d(x, y, bounds_error=False, fill_value='extrapolate')
        # Interpolate x and y values
        x_interp = np.linspace(min(x), max(x), 1000)
        y_interp = interp_func(x_interp)

        # Set inner boundary
        x1 = [1.5, 1.6, 1.8]
        y1 = [23, 21.7, 20.3]
        interp_func1 = interp1d(x1, y1, bounds_error=False, fill_value='extrapolate')
        x1_interp = np.linspace(min(x1), max(x1), 1000)
        y1_interp = interp_func1(x1_interp)

        x2 = [7, 6, 5, 4, 2.8, 1.8]
        y2 = [21.1, 20.7, 20.25, 20.15, 20.1, 20.3]
        interp_func2 = interp1d(x2, y2, bounds_error=False, fill_value='extrapolate')
        x2_interp = np.linspace(min(x2), max(x2), 1000)
        y2_interp = interp_func2(x2_interp)

        # Exclude values below the interpolated line
        excluded_df = df[(df[f814w_col] > interp_func(df[f475w_f814_col])) | (df[f814w_col] <= interp_func1(df[f475w_f814_col])) | (df[f814w_col] <= interp_func2(df[f475w_f814_col]))]
        ax0.scatter(excluded_df[f475w_f814_col], excluded_df[f814w_col], color='black')

        # Keep points within selection
        df = df[df[f814w_col] < interp_func(df[f475w_f814_col])]
        df = df[df[f814w_col] > interp_func1(df[f475w_f814_col])]
        df = df[df[f814w_col] > interp_func2(df[f475w_f814_col])]
        # ax0.scatter(excluded_df[f475w_f814_col], excluded_df[f814w_col], color='black')
        ax0.plot(x_interp, y_interp, color='black', linestyle='-', linewidth=1)
        ax0.plot(x1_interp, y1_interp, color='black', linestyle='-', linewidth=1)
        ax0.plot(x2_interp, y2_interp, color='black', linestyle='-', linewidth=1)

        # ax1
        scatter = ax1.scatter(df[f475w_f814_col], df[f814w_col], c=df[z_col], cmap='YlOrBr', s=1)
        ax1.plot(x_interp, y_interp, color='black', linestyle='-', linewidth=1)
        ax1.plot(x1_interp, y1_interp, color='black', linestyle='-', linewidth=1)
        ax1.plot(x2_interp, y2_interp, color='black', linestyle='-', linewidth=1)

        if save_rgb_selection == True:
            df.to_csv(save_filepath, index=False)
            print('Saved RGB selection to {}'.format(save_filepath))
        if savefile == True:
            plt.savefig(savefig_filepath)    
            print('Saved CMD to {}'.format(savefig_filepath))

    else:
        # Create the scatter plot
        scatter = plt.scatter(df[f475w_f814_col], df[f814w_col], c=df[z_col], cmap='YlOrBr', s=1)

        # Set the plot title and labels
        # plt.title('Scatter Plot')
        plt.xlabel(f475w_f814_col)
        plt.ylabel(f814w_col)

        # Invert the y-axis
        plt.gca().invert_yaxis()

        # Add colorbar
        cbar = plt.colorbar(scatter)
        cbar.set_label('MH')

        # Save the plot
        if savefile == True:
            plt.savefig(savefig_filepath)
                        
def photometry_spatial_map():
    """
    """

def photometry_binned_spatial_map():
    """
    """