import glob
import pandas as pd
import seaborn as sns

isochrone_path_list=glob.glob('/Users/mmckay/phd_projects/analysis_routine/DATA/CMD37_csvs/complete_tables/*.csv')
dfs = []
rgb_dfs = []

for isochrone_path in isochrone_path_list:
    df = pd.read_csv(isochrone_path)
    # Testing reduction of the ioschrone to isolate RGB stars
    df = df.loc[(df['F814W_appmag'] <= 23.0)]
    dfs.append(df)
    
    # Isolate RGB stars using the mean of the initial mass distribution
    rgb_df = df.loc[(df['Mini'] < df['Mini'].mean())]
    # Further reduction
    
    rgb_dfs.append(rgb_df)

    # save RGB files
    rgb_df.to_csv('/Users/mmckay/phd_projects/analysis_routine/DATA/CMD37_csvs/RGB_isochrone_tables/RGB_Zini_{}.csv'.format(df['Zini'].iloc[0]), index=False)
    

    

RGB_merged_df = pd.concat(rgb_dfs)
# RGB_merged_df = RGB_merged_df.loc[(RGB_merged_df['F475W_appmag-F814W_appmag'] <= 3.4)]
RGB_merged_df.to_csv('/Users/mmckay/phd_projects/analysis_routine/DATA/CMD37_csvs/RGB_isochrone_table.csv', index=False)


