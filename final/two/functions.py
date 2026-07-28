import numpy as np
import matplotlib.pyplot as plt
from gwosc.datasets import event_gps
import bilby
from bilby.core.prior import Uniform, PowerLaw, Cosine, Sine
from bilby.gw.conversion import convert_to_lal_binary_black_hole_parameters, generate_all_bbh_parameters
bilby.core.utils.log.setup_logger(log_level='WARNING')
from gwpy.timeseries import TimeSeries

def getData(name, interferrometer, duration):
    # setting variables
    #################################################################################################
    time = event_gps(name)
    start = time - (duration/2)
    end = time + (duration/2)

    sampling_frequency = 4096

    start_psd = time - 128
    end_psd = time + 128

    tempInterferrometer = bilby.gw.detector.get_empty_interferometer(interferrometer)
    
    # getting strain data
    #################################################################################################
    data = TimeSeries.fetch_open_data(
        interferrometer, start, end, sampling_frequency, cache=True)
    
    tempInterferrometer.set_strain_data_from_gwpy_timeseries(data)
    tempInterferrometer.sampling_frequency = sampling_frequency
    tempInterferrometer.duration = duration

    data.plot()
    plt.savefig(name + " " + interferrometer + " waveform.png", dpi = 150, bbox_inches='tight')

    # getting psd data
    #################################################################################################
    dataPsd = TimeSeries.fetch_open_data(
    interferrometer, start_psd, end_psd, sampling_frequency, cache=True)

    alpha = 2 * tempInterferrometer.strain_data.roll_off / duration
    psd = dataPsd.psd(fftlength = duration, overlap=0, window = ("tukey", alpha), method="median")

    tempInterferrometer.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(
    frequency_array= psd.frequencies.value, psd_array = psd.value)

    #graphing psd
    ##################################################################################################
    fig, ax = plt.subplots(figsize = [12,6])

    mask_psd = tempInterferrometer.power_spectral_density.frequency_array <= tempInterferrometer.maximum_frequency
    ax.loglog(tempInterferrometer.power_spectral_density.frequency_array[mask_psd],
        tempInterferrometer.power_spectral_density.asd_array[mask_psd])

    mask_strain = tempInterferrometer.strain_data.frequency_mask  
    ax.loglog(tempInterferrometer.strain_data.frequency_array[mask_strain],
        np.abs(tempInterferrometer.strain_data.frequency_domain_strain[mask_strain]))

    ax.set_xlim(10,1000)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel(r"Strain [strain/$\sqrt{Hz}$]")

    plt.savefig(name + " " + interferrometer + " psd waveform.png", dpi = 150, bbox_inches='tight')
    
    return tempInterferrometer

def originalPriors(time,prior):
    prior['chirp_mass'] = Uniform(name='chirp_mass', minimum=30.0,maximum=32.5)
    prior['mass_ratio'] = Uniform(name='mass_ratio', minimum=0.5, maximum=1)

    prior['a_1'] =  0.0
    prior['a_2'] =  0.0
    prior['tilt_1'] =  0.0
    prior['tilt_2'] =  0.0
    prior['phi_12'] =  0.0
    prior['phi_jl'] =  0.0

    prior['dec'] = Cosine(name="dec", minimum=-np.pi/2, maximum=np.pi/2)
    prior['ra'] = Uniform(name="ra", minimum=0, maximum=2*np.pi, boundary='periodic')
    prior['theta_jn'] =  1.89694 #fixed to the one in the tutorial?
    prior['psi'] =  0.532268 #fixed to the one in the tutorial?
    prior['phase'] = Uniform(name="phase", minimum=0, maximum=2*np.pi)
    prior['geocent_time'] = Uniform(name="geocent_time", minimum=time-0.1, maximum=time+0.1)
    prior['luminosity_distance'] = PowerLaw(alpha=2, name='luminosity_distance', minimum=50, maximum=2000, unit='Mpc', latex_label='$d_L$')
    return prior

def maximalPriors(time,prior, minimum, maximum):
    prior['chirp_mass'] = Uniform(name='chirp_mass', minimum=minimum, maximum=maximum)
    prior['mass_ratio'] = Uniform(name='mass_ratio', minimum=0.5, maximum=1)
    
    prior['a_1'] =  Uniform(name = 'a_1', minimum = 0, maximum = 0.8)
    prior['a_2'] =  Uniform(name = 'a_2', minimum = 0, maximum = 0.8)
    prior['tilt_1'] =  Sine(name = 'tilt_1', minimum = 0, maximum = np.pi)
    prior['tilt_2'] =  Sine(name = 'tilt_2', minimum = 0, maximum = np.pi)
    prior['phi_12'] =  Uniform(name = 'phi_12', minimum = 0, maximum = 2 * np.pi)
    prior['phi_jl'] =  Uniform(name = 'phi_jl', minimum = 0, maximum = 2 * np.pi)
    prior['dec'] = Cosine(name="dec", minimum=-np.pi/2, maximum=np.pi/2)
    prior['ra'] = Uniform(name="ra", minimum=0, maximum=2*np.pi, boundary='periodic')
    prior['theta_jn'] =  Sine(name = 'theta_jn', minimum = 0, maximum = np.pi)
    prior['psi'] =  Uniform(name = 'psi', minimum = 0, maximum = np.pi)
    
    prior['phase'] = Uniform(name="phase", minimum=0, maximum=2*np.pi)
    prior['geocent_time'] = Uniform(name="geocent_time", minimum=time-0.1, maximum=time+0.1)
    prior['luminosity_distance'] = PowerLaw(alpha=2, name='luminosity_distance', minimum=50, maximum=2000, unit='Mpc', latex_label='$d_L$')
    return prior

def experimentalPriorsPartialSpin(time, prior,minimum, maximum):
    #MASSES
    prior['chirp_mass'] = Uniform(name='chirp_mass', minimum=minimum, maximum=maximum)
    prior['mass_ratio'] = Uniform(name='mass_ratio', minimum=0.5, maximum=1)

    #SPINS
    prior['a_1'] =  Uniform(name = 'a_1', minimum = 0, maximum = 0.8)
    prior['a_2'] =  Uniform(name = 'a_2', minimum = 0, maximum = 0.8)
    prior['tilt_1'] = 0.0
    prior['tilt_2'] = 0.0
    prior['phi_12'] = 0.0
    prior['phi_jl'] = 0.0

    prior['ra'] = 0.0  #these make it SUPER long...
    prior['dec'] = 0.0  
    prior['theta_jn'] = Sine(name='theta_jn', minimum=0, maximum=np.pi)
    prior['psi'] =  Uniform(name = 'psi', minimum = 0, maximum = np.pi)
    prior['phase'] = Uniform(name="phase", minimum=0, maximum=2*np.pi) 
    prior['luminosity_distance'] = PowerLaw(alpha=2, name='luminosity_distance', minimum=50, maximum=5000, unit='Mpc', latex_label='$d_L$')
    prior['geocent_time'] = Uniform(name="geocent_time", minimum=time-0.1, maximum=time+0.1)

    return prior

def experimentalPriorNoSpin(time, prior, minimum, maximum):
    prior['chirp_mass'] = Uniform(name='chirp_mass', minimum=minimum, maximum=maximum)
    prior['mass_ratio'] = Uniform(name='mass_ratio', minimum=0.5, maximum=1)
    
    prior['a_1'] = 0.0 # 0 to 1
    prior['a_2'] = 0.0
    prior['tilt_1'] = 0.0 # cause precession,set to 0 for aligns spins
    prior['tilt_2'] = 0.0 
    prior['phi_12'] = 0.0 
    prior['phi_jl'] = 0.0
    
    prior['dec'] = Cosine(name="dec", minimum=-np.pi/2, maximum=np.pi/2)
    prior['ra'] = Uniform(name="ra", minimum=0, maximum=2*np.pi, boundary='periodic')
    prior['theta_jn'] = Sine(name='theta_jn', minimum=0, maximum=np.pi) 
    prior['psi'] = Uniform(name='psi', minimum=0, maximum=np.pi, boundary='periodic')
    prior['phase'] = Uniform(name="phase", minimum=0, maximum=2*np.pi)
    prior['geocent_time'] = Uniform(name="geocent_time", minimum=time-0.1, maximum=time+0.1)
    prior['luminosity_distance'] = PowerLaw(alpha=2, name='luminosity_distance', minimum=50, maximum=5000, unit='Mpc')
    return prior

def experimentalPriorSpin(time, prior, minimum, maximum):
    prior['chirp_mass'] = Uniform(name='chirp_mass', minimum=minimum, maximum=maximum)
    prior['mass_ratio'] = Uniform(name='mass_ratio', minimum=0.5, maximum=1)
    
    prior['a_1'] =  Uniform(name = 'a_1', minimum = 0, maximum = 0.8)
    prior['a_2'] =  Uniform(name = 'a_2', minimum = 0, maximum = 0.8)
    prior['tilt_1'] = 0.0 # cause precession, set to 0 for aligns spins hahahah still probably wrong but whatever
    prior['tilt_2'] = 0.0 
    prior['phi_12'] = 0.0 
    prior['phi_jl'] = 0.0
    
    prior['dec'] = Cosine(name="dec", minimum=-np.pi/2, maximum=np.pi/2)
    prior['ra'] = Uniform(name="ra", minimum=0, maximum=2*np.pi, boundary='periodic')
    prior['theta_jn'] = Sine(name='theta_jn', minimum=0, maximum=np.pi) 
    prior['psi'] = Uniform(name='psi', minimum=0, maximum=np.pi, boundary='periodic')
    prior['phase'] = Uniform(name="phase", minimum=0, maximum=2*np.pi)
    prior['geocent_time'] = Uniform(name="geocent_time", minimum=time-0.1, maximum=time+0.1)
    prior['luminosity_distance'] = PowerLaw(alpha=2, name='luminosity_distance', minimum=50, maximum=5000, unit='Mpc')
    return prior

def sampler(prior, nlive, dlogz, detectors, name):
    
    waveform_arguments = dict(
        waveform_approximant='IMRPhenomXP', reference_frequency=100., catch_waveform_errors=True)

    waveform_generator = bilby.gw.WaveformGenerator(
        frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
        waveform_arguments=waveform_arguments,
        parameter_conversion=convert_to_lal_binary_black_hole_parameters)

    likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
        detectors, waveform_generator, priors=prior,
        time_marginalization=True, phase_marginalization=True, distance_marginalization=True)

    result_short = bilby.run_sampler(
        likelihood, prior, sampler='dynesty', outdir='short', label=name,
        conversion_function=bilby.gw.conversion.generate_all_bbh_parameters,
        nlive=nlive, 
        dlogz= dlogz,  
        clean=True,)
    
    return result_short

