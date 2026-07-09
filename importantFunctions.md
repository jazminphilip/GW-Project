# from gwosc.datasets
find_datasets(type) just outputs a catalog based on what is the type 
run_segment
event_gps(name) gives gps time

# from gwpy.timeseries import TimeSeries
TimeSeries.fetch_open_data(interferometer, start gps time, end gps time, cache = True) just downloads data onto a variable
specgram = ldata.spectrogram2(fftlength=4, overlap=2, window='hann') ** (1/2.) TO GET Q TRANFORM
plot1psd = ldata1.psd(fftlength=4, overlap=2, window='hann') TO GET PSD

# from pycbc.waveform import get_td_waveform
hp, hc = figures = get_td_waveform(approximant="IMRPhenomT", mass1 = x , mass2 = y, inclination=0, delta_t=1/2048, f_lower = 30) this is to make an ideal wave form
estimated_psd = pycbc.psd.welch(ts,seg_len=seg_len,seg_stride=seg_stride)


# import matplotlib.pyplot as plt
plt.figure(figsize=plt.figaspect(0.4))
plt.plot(hp.sample_times, hp, label='Plus Polarization')
plt.plot(hp.sample_times, hc, label='Cross Polarization')\
plt.legend()
plt.grid()

fig, ax = plt.subplots(1,1, figsize=[6,4])
ax.loglog(freqs_welch, psd_welch, label='Welch est.')
ax.legend()

plt.hist(data_whitened, bins=100, density=True, alpha=0.7, color='lightcoral', edgecolor='black')


# import scipy
freqs_welch, psd_welch = scipy.signal.welch(ldata.value, fs,  window='hann', nperseg=4096)

# import numpy as np
cross_correlation = numpy.zeros([len(data)-len(hp1)]) makes zero array
hp1_numpy = hp1.numpy() converts to numpy array
data_whitened = (ts.to_frequencyseries() / psd**0.5).to_timeseries()
snr_series = np.correlate(datan, hp1n, mode='same') datan and cross correlate with a function

# from pycbc.psd import interpolate, inverse_spectrum_truncation
# from pycbc.filter import matched_filter
psd = interpolate(psd, conditioned.delta_f) to make the x values of psd and delta_f match
psd = inverse_spectrum_truncation(psd, int(4 * conditioned.sample_rate),
                                  low_frequency_cutoff=15)
snr = matched_filter(template, conditioned,
                     psd=psd, low_frequency_cutoff=20)

