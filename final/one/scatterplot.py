# you have to know csvfile manipulation through panda
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("events.csv")
df = df.dropna(subset=["mass_1_source", "mass_2_source", "chirp_mass_source"])

name = "GW150914"

event = df[df["name"].str.contains(name)]

m1 = event["mass_1_source"].values[0] # literally the name on the csv file
m2 = event["mass_2_source"].values[0]

fig, ax = plt.subplots(figsize=(7, 6))

ax.axline(xy1=(0,0), slope=1.0, color='black', linestyle='--')

ax.scatter(
    df["mass_1_source"], #lists an array and graphs it
    df["mass_2_source"],
    s=20, # point size
    color="lightgray",
    edgecolor="gray",
    linewidth=0.3,
    label="GWTC catalog events",
    zorder=2,
)

ax.scatter(
    m1,
    m2,
    s=200,
    color="crimson",
    edgecolor="black",
    linewidth=1.2,
    marker="*",
    label= name + "\nM1 = " + str(m1) + ", M2 = " + str(m2),
    zorder=3,
)

max_mass = df["mass_1_source"].max()

ax.set_xlabel("Primary mass, $m_1$ ($M_\\odot$)", fontsize=12)
ax.set_ylabel("Secondary mass, $m_2$ ($M_\\odot$)", fontsize=12)
ax.set_title("Component Masses Across the GWTC Catalog", fontsize=13)
ax.legend(frameon=False, fontsize=10)
ax.set_xlim(0, max_mass)
ax.set_ylim(0, max_mass)

plt.savefig("scatterplot.png", dpi=150)

