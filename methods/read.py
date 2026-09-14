import pandas as pd
import corner
import matplotlib.pyplot as plt

df = pd.read_csv('GW241109_115924 posterior samples.csv')

# Select only the 5 columns you want
columns = ['mass_1','mass_2', 'chirp_mass', 'mass_ratio']  # change these to what you want
df_subset = df[columns]

# Convert to numpy and plot
data = df_subset.to_numpy()
fig = corner.corner(data, 
                    labels=columns,
                    show_titles=True,
                    title_fmt='.3f',
                    quantiles=[0.05, 0.5, 0.95])  # 90% CI

plt.show()
plt.savefig("GW241109_115924 corners from csv.png", dpi = 500)