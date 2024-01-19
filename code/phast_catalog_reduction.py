import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px

df = pd.read_csv('/Users/mmckay/phd_projects/analysis_routine/DATA/phast_f475W_f814W_table.csv')

# Gregersen et al 2015 criteria for RGB stars
# (df['f475w_flag'] == 0.0) & (df['f814w_flag'] == 0.0) DOLPHOT flags... flags < 8 are good
phast_df = df.loc[(df['f475w_snr'] >= 4.0) & (df['f814w_snr'] >= 4.0)
                     & (df['f475w_crowd']+df['f814w_crowd'] <= 1.0)
                     & ((df['f475w_sharp']+df['f814w_sharp'])**2 <= 0.075)
                    #  & (df['f814w_vega'] <= 23.0) # Apply after foreground reddening correction
                     ].copy()

# Apply foreground reddening correction
# A_lam/A_V: F814W=0.59696, F475W=1.21182 from CMD 3.7 output A_V = 0.17
phast_df.loc[:, 'f814w_vega_ecorr'] = phast_df['f814w_vega'] + 0.596*0.17
phast_df.loc[:, 'f475w_vega_ecorr'] = phast_df['f475w_vega'] + 1.212*0.17

# Calculate F475W-F814W
phast_df.loc[:, 'f475w-f814w_ecorr'] = phast_df['f475w_vega_ecorr']-phast_df['f814w_vega_ecorr']
phast_df.loc[:, 'f475w-f814w_vega'] = phast_df['f475w_vega']-phast_df['f814w_vega']


# Apply brightness cut to the entire catalog after foreground reddening correction
phast_df = phast_df.loc[(phast_df['f814w_vega_ecorr'] <= 23.0)]

# Save the dataframe as a CSV file
print('Saving the dataframe as a CSV file...')
phast_df.to_csv('/Users/mmckay/phd_projects/analysis_routine/DATA/reduced_phast_catalog.csv', index=False)

# # Interactive pairplot with hover tool
# fig = px.scatter_matrix(phast_df, dimensions=['f475w_vega_ecorr', 'f814w_vega_ecorr', 'f475w-f814w_ecorr', 'f475w_snr', 'f814w_snr'],
#                         hover_data=phast_df.columns)
# fig.show()

# fig.write_html("/Users/mmckay/phd_projects/analysis_routine/FIGURES/phast_scatter_matrix_catalog.html")

# # Save the figure as JPEG
# fig.write_html("/Users/mmckay/phd_projects/analysis_routine/FIGURES/phast_scatter_matrix_catalog.jpg", format='jpeg')


