# pre_processing scripts for processing CMD isochrone tables

import pandas as pd
import os
import re
from tqdm import tqdm
import glob
import numpy as np

def read_sort_isochrone_table(filename, output_dir):
    
    """
    Reads in the CMD table line by line and sorts the isochrone interval into separate dataframes and 
    saves them as csvs to an output directory

    Parameters
    ----------
    filename : str
        Path to the isochrone the complete isochrone table produced by the CMD3.7 isochrone generator (http://stev.oapd.inaf.it/cgi-bin/cmd)
    output_dir : str
        Path to the directory where the output csvs will be saved(i.e /Users/mmckay/phd_projects/analysis_routine/DATA/CMD37_csvs/)

    """
    data_frames = []
    # Initialize variables
    table_start_pattern = (
        r"# Zini     MH  (.+)"  # Regular expression pattern to identify table start
    )

    # Open the file for reading
    table_num = 0
    with open(filename, "r") as file:
        for line in tqdm(file, desc="Processing file", unit="line"):
            # line = file.readline()
            # print(f'Line: {line}')
            # Looks for line that matches the pattern for every new column
            match = re.match(table_start_pattern, line)
            if match:
                # If the column matches the pattern, then start a new table
                table_num += 1
                stripped_line = line.strip()
                stripped_line = stripped_line.replace("\t", " ")
                words = stripped_line.split()
                space_separated_line = " ".join(words[1:])
                comma_separated_line = space_separated_line.replace(" ", ",")
                list_of_colname = comma_separated_line.split(",")

                
                if table_num > len(data_frames):
                    # print(f"{table_num} New dataframe created")
                    df = pd.DataFrame(columns=list_of_colname)
                    data_frames.append(df)
                else:
                    df = data_frames[table_num]
            # If first element on the line is a comment then ignore it
            elif line[0] == "#":
                pass
            # if line doesnt meet any of the above criteria then it is assumed to be isochrone row data and is added to the newest dataframe
            else:
                stripped_line = line.strip()
                stripped_line = stripped_line.replace("\t", " ")
                words = stripped_line.split()
                space_separated_line = " ".join(words)
                comma_separated_line = space_separated_line.replace(" ", ",")
                list_of_data = comma_separated_line.split(",")
                list_of_data = [float(i) for i in list_of_data]

                # Add data to dataframe
                df.loc[len(df)] = list_of_data



    # Save each isochrone as a separate csv
    os.makedirs(
        f"{output_dir}/{filename.split('/')[-1].split('.')[0]}",
        exist_ok=True,
    )


    ##! Issue: Some dataframes are being merged with next isochrone randomly 
    # for i, df in enumerate(data_frames):
    #     df.to_csv(f"{output_dir}/{filename.split('/')[-1].split('.')[0]}/{i}_{filename.split('/')[-1].split('.')[0]}.csv", index=False)

    # Merge isochrones together to single csv file
    merged_df = pd.concat(data_frames)
    merged_df.to_csv(
        f"{output_dir}/{filename.split('/')[-1].split('.')[0]}/merged_{filename.split('/')[-1].split('.')[0]}.csv",
        index=False,
    )

    print("Merged Dataframe Complete")
    # Split the newly merged by metallicity value and replace saved dataframe
    # Split cmd_df by metallicity value
    cmd_df_split = pd.DataFrame()
    for zini_value in merged_df['Zini'].unique():
        cmd_df_split = pd.DataFrame(merged_df[merged_df['Zini'] == zini_value])
        cmd_df_split.to_csv('{}/{}/output_{}.csv'.format(output_dir, filename.split('/')[-1].split('.')[0], zini_value), index=False)

    return merged_df

def cmd_isochrone_processing(isochrone_csv, distance_modules=24.45, save_tables=True):
    """
    Processes the isochrone csv produced by read_sort_isochrone_table by
    - Converting photmetry to from absoluete to apprent magnitude using the distance modules
    - calulate the difference of F475W and F814 for color magnitude diagram for both magnituides
    - Applying Foreground color extinction correction

    Parameters
    ----------
    isochrone_csv : str
        Path to the isochrone csv produced by read_sort_isochrone_table
    distance_modules : float, optional
        Distance modules to the Androneda Galaxy, by default 24.45 (Gregersien et al. 2019)
    save_tables : bool, optional
    
    """

    # Convert absolute magnitude to apparent magnitude
    cmd_df = pd.read_csv(isochrone_csv)
    d = 10 ** ((distance_modules + 5) / 5)
    cmd_df['F475W_appmag'] = cmd_df['F475Wmag'] - 5 + 5 * np.log10(d)
    cmd_df['F814W_appmag'] = cmd_df['F814Wmag'] - 5 + 5 * np.log10(d)
    cmd_df['F475Wmag-F814Wmag'] = cmd_df['F475Wmag'] - cmd_df['F814Wmag']
    cmd_df['F475W_appmag-F814W_appmag'] = cmd_df['F475W_appmag'] - cmd_df['F814W_appmag']

    # save the processed isochrone table
    if save_tables is True:
        cmd_df.to_csv(f'{isochrone_csv}')

    return cmd_df



def isochrone_rgb_selection(isochrone_dirpath, output_dir):
    """
    Selects RGB stars from an isochrone table by 2 criteria
    - F814W_appmag <= 23.0
    - inital mass <= mean of the intial mass of AGB and RGB stars
    
    Saves them as individual csv for each metallicity interval and then a single complete file
    
    Parameters
    ----------
    isochrone_path : str
        Path to the isochrone directory after it has been processed by read_sort_isochrone_table where
        each metallicity interval is a separate csv(i.e. output)
    output_dir : str
        Path to the directory where the output csvs will be saved
        (i.e /Users/mmckay/phd_projects/analysis_routine/DATA/CMD37_csvs/complete_tables/)

    Returns
    -------
    RGB_merged_df : pandas dataframe
        Dataframe of all RGB stars from all metallicity intervals from isochrone table

    """

    isochrone_path_list=glob.glob(f'{isochrone_dirpath}/*/output*.csv')
    dfs = []
    rgb_dfs = []

    for isochrone_path in isochrone_path_list:
        df = cmd_isochrone_processing(isochrone_path)
        # df = pd.read_csv(isochrone_path)
        # Testing reduction of the ioschrone to isolate RGB stars
        df = df.loc[(df['F814W_appmag'] <= 23.0)]
        dfs.append(df)
        
        # Isolate RGB stars using the mean of the initial mass distribution
        rgb_df = df.loc[(df['Mini'] < df['Mini'].mean())]
        # Further reduction
        
        rgb_dfs.append(rgb_df)

        # save RGB files
        os.makedirs(f'{output_dir}/RGB/', exist_ok=True)
        os.makedirs(f'{output_dir}/RGB_AGB/', exist_ok=True)
        rgb_df.to_csv('{}/RGB/RGB_Zini_{}.csv'.format(output_dir, df['Zini'].iloc[0]), index=False)
        df.to_csv('{}/RGB_AGB/RGB_AGB_Zini_{}.csv'.format(output_dir, df['Zini'].iloc[0]), index=False)
        
    RGB_merged_df = pd.concat(rgb_dfs)
    RGB_AGB_merged_df = pd.concat(dfs)

    #  Save dataframes as csvs
    RGB_merged_df.to_csv(f'{output_dir}/RGB/RGB_isochrone_table.csv', index=False) 
    RGB_AGB_merged_df.to_csv(f'{output_dir}/RGB_AGB/RGB_AGB_isochrone_table.csv', index=False)

    return RGB_merged_df
