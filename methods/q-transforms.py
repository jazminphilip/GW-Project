import gwpy
from gwosc.datasets import event_gps
from gwpy.timeseries import TimeSeries
import matplotlib.pyplot as plt

name = "GW250120_042414"
interferrometers = {'L1'}
    #'H1','L1','V1','G1','K1'}

gps = event_gps(name)
lower_bound = int(gps) - 512
upper_bound = int(gps) + 512

for interferrometer in interferrometers:
    try:    
        data = TimeSeries.fetch_open_data(interferrometer, lower_bound, upper_bound, verbose = True, cache = True)
        freq_range = (20, 1000)
        q_range = (16, 24)

        gated_data = data.gate(tzero = 0.25, tpad = 0.25)
        q_transform_gated_data = gated_data.q_transform(frange = freq_range, qrange = q_range, outseg=(gps-1, gps+1)) # OUTSEG IS LIKE WHAT IS FOCUSED ON
        plot = q_transform_gated_data.plot()

        ax = plot.gca()
        ax.set_epoch(gps)

        ax.set_yscale('log')
        ax.set_xlim(gps-.5,gps+.5)

        ax.colorbar(label="Normalised energy")
        ax.colorbar().mappable.set_clim(0,20)
        ax.set_title(name + " " + interferrometer + " Q Transform")
        plt.savefig(name + " " + interferrometer + " Q Transform.png", dpi=150, bbox_inches='tight')   
    except:
        print(interferrometer + " not available.") 