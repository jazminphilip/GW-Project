import numpy as np
import matplotlib.pyplot as plt
from gwpy.timeseries import TimeSeries
from pycbc.waveform import get_td_waveform
from pycbc.filter import matched_filter
from pycbc.psd import interpolate, inverse_spectrum_truncation
import pycbc.types
from gwosc.datasets import event_gps

# ── 1. Fetch real strain data from both detectors ──────────────────────
name = "GW231206_010629"
event_gps = event_gps(name)
detectors = ["H1", "L1"]
f_lower = 30.0

data = {}
psds = {}

for det in detectors:
    print(f"Fetching {det} strain data around " + name + "...")
    strain = TimeSeries.fetch_open_data(
        det, event_gps - 32, event_gps + 32, sample_rate=4096, cache=True
    )
    ts = pycbc.types.TimeSeries(strain.value, delta_t=strain.dt.value, epoch=strain.t0.value)
    data[det] = ts

    seg_len_psd = 4
    psd = ts.psd(seg_len_psd)
    psd = interpolate(psd, ts.delta_f)
    psd = inverse_spectrum_truncation(
        psd, int(seg_len_psd * ts.sample_rate), low_frequency_cutoff=f_lower
    )
    psds[det] = psd

# ── 2. Set up the finer 2D grid of trial (m1, m2) templates ────────────
m1_grid = np.linspace(15, 50, 60)
m2_grid = np.linspace(15, 45, 60)

snr_grid = np.full((len(m1_grid), len(m2_grid)), np.nan)

print("Running 2D matched-filter grid search across H1 + L1 "
      f"({len(m1_grid)}x{len(m2_grid)} = {len(m1_grid)*len(m2_grid)} templates)...")
print("This will take longer than the single-detector version -- be patient.")

for i, m1 in enumerate(m1_grid):
    for j, m2 in enumerate(m2_grid):
        try:
            combined_snr_sq = 0.0
            valid = True
            for det in detectors:
                hp, _ = get_td_waveform(
                    approximant="IMRPhenomD",
                    mass1=m1, mass2=m2,
                    delta_t=data[det].delta_t, f_lower=f_lower,
                )
                if len(hp) > len(data[det]):
                    valid = False
                    break
                hp.resize(len(data[det]))
                hp = hp.cyclic_time_shift(hp.start_time)

                snr = matched_filter(
                    hp, data[det], psd=psds[det], low_frequency_cutoff=f_lower
                )
                crop = int(4 * data[det].sample_rate)
                snr = snr[crop:-crop]
                peak_snr = abs(snr).numpy().max()

                # Combine detectors in quadrature -- this is the network SNR
                combined_snr_sq += peak_snr ** 2

            if valid:
                snr_grid[i, j] = np.sqrt(combined_snr_sq)
        except Exception:
            snr_grid[i, j] = np.nan
    if i % 10 == 0:
        print(f"  m1 = {m1:.1f} done ({i+1}/{len(m1_grid)} rows)")

print("Grid search complete.")

# ── 3. Identify the best-fit point ──────────────────────────────────────
best_idx = np.unravel_index(np.nanargmax(snr_grid), snr_grid.shape)
best_m1 = m1_grid[best_idx[0]]
best_m2 = m2_grid[best_idx[1]]
best_snr = snr_grid[best_idx]

best_chirp_mass = (best_m1 * best_m2) ** (3/5) / (best_m1 + best_m2) ** (1/5)

print(name)
print(f"\nBest-fit grid point: m1 = {best_m1:.2f}, m2 = {best_m2:.2f} solar masses")
print(f"Network SNR at best-fit point: {best_snr:.1f}")
print(f"Sampled Chirp Mass: {best_chirp_mass:.2f} solar masses")

# ── 4. Plot the SNR heatmap ─────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 6))
mesh = ax.pcolormesh(m2_grid, m1_grid, snr_grid, shading="auto", cmap="viridis")
plt.colorbar(mesh, label="Network matched-filter SNR (H1 + L1)")

ax.scatter(
    best_m2, best_m1, color="red", marker="*", s=250,
    edgecolor="white", linewidth=1.2, label="Grid best-fit point",
)
# ax.scatter(
#    29.0, 36.0, color="white", marker="x", s=120,
#    linewidth=2.5, label="Published GWTC value",
#)

ax.set_xlabel("$m_2$ ($M_\\odot$)", fontsize=12)
ax.set_ylabel("$m_1$ ($M_\\odot$)", fontsize=12)
ax.set_title("2D Matched-Filter Network SNR Grid — " + name + " (real data)", fontsize=12)
ax.legend(fontsize=10)

plt.tight_layout()
plt.savefig(name + " grid search.png", dpi=150)
plt.show()