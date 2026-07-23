# you have to know csvfile manipulation through panda
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("events.csv")
df = df.dropna(subset=["mass_1_source", "mass_2_source", "chirp_mass_source"])

name = "GW150914"

event = df[df["name"].str.contains(name)]

chirp_mass = event["chirp_mass_source"].item()

fig, ax = plt.subplots(figsize=(7, 5))

ax.hist(
    df["chirp_mass_source"],
    bins=30,
    color="lightgray",
    edgecolor="gray",
    linewidth=0.5,
    label="GWTC catalog events",
)

ax.axvline(
    chirp_mass,
    color="crimson",
    linewidth=2.5,
    label = name +  " (" + str(chirp_mass) + " $M_\\odot$)",
)

ax.set_xlabel("Chirp mass, $\\mathcal{M}$ ($M_\\odot$)", fontsize=12)
ax.set_ylabel("Number of events", fontsize=12)
ax.set_title("Chirp Mass Distribution Across the GWTC Catalog", fontsize=13)
ax.legend(frameon=False, fontsize=10)

plt.tight_layout()
plt.savefig("histogram.png", dpi=150)
plt.show()
