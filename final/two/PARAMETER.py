import numpy as np
import matplotlib.pyplot as plt
from gwosc.datasets import event_gps
import bilby
from bilby.core.prior import Uniform, PowerLaw, Cosine
from bilby.gw.conversion import convert_to_lal_binary_black_hole_parameters, generate_all_bbh_parameters

bilby.core.utils.log.setup_logger(log_level='WARNING')

from gwpy.timeseries import TimeSeries

import functions

name = "GW231206_010629"
time = event_gps(name)
duration = 8

# interferrometers = ['G1','H1', 'L1', 'V1', 'K1']

interferrometers = ['H1', 'L1']
detectors = []

for det in interferrometers:
    detectors.append(functions.getData(interferrometer= det, name = name, duration = duration))

#################################################################

prior = bilby.core.prior.PriorDict()
functions.experimentalPriorNoSpin(time = time ,prior = prior, minimum = 10, maximum = 20)

##################################################################

result_short = functions.sampler(prior = prior, nlive = 5000, dlogz = .2, detectors = detectors, name = name)

##################################################################
result_short.posterior
result_short.posterior.to_csv(name + " posterior samples.csv", index=False)

###################################################################
# parametersfull = ['chirp_mass','mass_1','mass_2','a_1','a_2','tilt_1','tilt_2','phi_12','phi_jl','dec','ra','theta_jn','psi','phase','geocent_time','luminosity_distance']

parametersfull = ['mass_1','mass_2', 'chirp_mass', 'mass_ratio']
dictionary = {}

parameters = dict()
for x in parametersfull:
    if np.quantile(result_short.posterior[x], .5) != 0.0:    
        y = float(np.quantile(result_short.posterior[x].values, .5))
        if not (x == 'a_1' or x == 'a_1'):
            parameters[x] = y
    print(x + ": " + str(y))
fig = result_short.plot_corner(parameters=parameters)
fig.savefig(name + " corner plots full.png", dpi=150, bbox_inches='tight')
