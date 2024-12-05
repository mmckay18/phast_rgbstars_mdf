from scipy.interpolate import LinearNDInterpolator
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np



RGB_merged_df = pd.read_csv('/Users/mmckay/phd_projects/analysis_routine/DATA/CMD37_csvs/RGB_isochrone_table.csv')
phast_df = pd.read_csv('/Users/mmckay/phd_projects/analysis_routine/DATA/reduced_phast_catalog.csv')


# Create a 2D array of F814W and F475W-F814W magnitudes from RGB_merged_df
cmd_arr = RGB_merged_df[['F814W_appmag', 'F475W_appmag-F814W_appmag']].values

# Mini,int_IMF,Mass,logL,logTe,logg,
# Create a LinearNDInterpolator object with 'MH' as the values
interp_MH = LinearNDInterpolator(cmd_arr, RGB_merged_df['MH'], fill_value=-99.0)
interp_Mini = LinearNDInterpolator(cmd_arr, RGB_merged_df['Mini'], fill_value=-99.0)
interp_int_IMF = LinearNDInterpolator(cmd_arr, RGB_merged_df['int_IMF'], fill_value=-99.0)
interp_Mass = LinearNDInterpolator(cmd_arr, RGB_merged_df['Mass'], fill_value=-99.0)
interp_logL = LinearNDInterpolator(cmd_arr, RGB_merged_df['logL'], fill_value=-99.0)
interp_logTe = LinearNDInterpolator(cmd_arr, RGB_merged_df['logTe'], fill_value=-99.0)
interp_logg = LinearNDInterpolator(cmd_arr, RGB_merged_df['logg'], fill_value=-99.0)



# Create a 2D array of f814w_vega and f475w-f814w magnitudes from phast_df
phast_df['f475w_vega-f814w_vega'] = phast_df['f475w_vega'] - phast_df['f814w_vega']

# PHAST catalog color magnitude diagram values
filtered_arr = phast_df[['f814w_vega', 'f475w_vega-f814w_vega']].values

# Interpolate the 'MH' values for each row in phast_df
interpolated_MH = interp_MH(filtered_arr)
interpolated_Mini = interp_Mini(filtered_arr)
interpolated_int_IMF = interp_int_IMF(filtered_arr)
interpolated_Mass = interp_Mass(filtered_arr)
interpolated_logL = interp_logL(filtered_arr)
interpolated_logTe = interp_logTe(filtered_arr)
interpolated_logg = interp_logg(filtered_arr)

# Add the interpolated 'MH' values as a new column in phast_df
phast_df['interpolated_MH'] = interpolated_MH
phast_df['interpolated_Mini'] = interpolated_Mini
phast_df['interpolated_int_IMF'] = interpolated_int_IMF
phast_df['interpolated_Mass'] = interpolated_Mass
phast_df['interpolated_logL'] = interpolated_logL
phast_df['interpolated_logTe'] = interpolated_logTe
phast_df['interpolated_logg'] = interpolated_logg

print(phast_df.shape)
print(phast_df.head(15))

phast_df = phast_df[phast_df['interpolated_MH'] != -99.0].dropna()
print(phast_df.shape)
print(phast_df.head(15))

# Drop  M/H values > 0.6 ~ CHECK IF THIS THE RIGHT THING TO DO 
# phast_df = phast_df[phast_df['interpolated_MH'] <= 0.6].dropna()

# Save the merged dataframe as CSV
phast_df.to_csv('/Users/mmckay/phd_projects/analysis_routine/DATA/interpolated_phast_MH_catalog.csv', index=False)


# add the closest index as a new column in phast_df
# phast_df['closest_index'] = closest_index.astype(int)


# merge phast_df with RGB_merged_df on closest_index
# RGB_metallicity_df = pd.merge(phast_df, RGB_merged_df, left_on='closest_index', right_index=True)
# RGB_metallicity_df = RGB_metallicity_df[RGB_metallicity_df['closest_index'] != 999999999]


# print the merged dataframe
# RGB_metallicity_df[['f475w_vega', 'f814w_vega', 'F475W_appmag', 'F814W_appmag', 'f475w-f814w_ecorr']].describe()

# Save the merged dataframe as CSV
# RGB_metallicity_df.to_csv('/Users/mmckay/phd_projects/analysis_routine/DATA/phast_RGB_MH_table.csv', index=False)

# Plot the metallicity distribution
# joint = sns.jointplot(x=phast_df['f814w_vega-f475w_vega'], y=phast_df['f814w_vega'], kind='scatter', color='grey', hue=phast_df['interpolated_MH'], height=0.5)
# plt.ylim(24, 19)
# joint.ax_marg_x.set_xlim(-1, 8)
# plt.savefig("/Users/mmckay/phd_projects/analysis_routine/FIGURES/phast_metallicity_distribution.jpg", format='jpeg')
# plt.close()

# sns.histplot(data=phast_df, x='interpolated_MH', label=f'N={len(phast_df)}', bins=20, kde=True, stat='count')
# plt.legend()
# plt.savefig("/Users/mmckay/phd_projects/analysis_routine/FIGURES/phast_mdf.jpg", format='jpeg')
# plt.close()

# sns.scatterplot(data=phast_df, x='ra', y='dec', hue='interpolated_MH', s=1, palette='magma', label=f'N={len(phast_df)}')
# plt.legend()
# plt.savefig("/Users/mmckay/phd_projects/analysis_routine/FIGURES/phast_mdf.jpg", format='jpeg')
# plt.close()

