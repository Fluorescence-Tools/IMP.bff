"""
Fluorescence decay pile-up
==========================
In single photon counting at high count rates subseqentially emitted 
photons cannot be detected. This effect is called pile up, as it changes
the the shape of fluorescence decay curves. This examples illustrates 
how pile-up changes the shape of fluorescence decays.
"""
import pylab as plt
import numpy as np
import scipy.stats
import IMP.bff

n_channels = 128
irf_position, irf_width = 1.0, 0.2
time_axis = np.linspace(0, 10, n_channels)
dt = time_axis[1] - time_axis[0]

irf = scipy.stats.norm.pdf(time_axis, loc=irf_position, scale=irf_width)
lifetime_spectrum = np.array([100., 4])
data = np.zeros_like(irf)
stop = len(irf)
start = 0

IMP.bff.decay_fconv(fit=data, irf=irf, x=lifetime_spectrum, start=start, stop=stop, dt=dt)


# %%
# The number the model is scaled to a number of photons typically recorded in a eTCSPC 
# experiment. In an eTCSPC experiment usually 1e6 - 20e6 photons are recorded. 
#
# The pile up depends on the instrument dead time, the repetition rate used to excite
# the sample and the total number of registered photons. Here, to highlight the effect
# of pile up an unrealistic combination of the measurement time and the number of photons 
# is used. In modern counting electronics the dead time is often around 100 nano
# seconds.
# 
# In this example there is no data. Thus, the model without pile up is used as
# "data". In a real experiment use the experimental histogram in the data argument.

period = 8.0
pile_up_parameter = {
    'repetition_rate': 1./period * 1000.,  # Repetition rate is in MHz
    'instrument_dead_time': 120.,  # Instrument dead time nano seconds
    'measurement_time': 0.001,  # Measurement time in seconds
    'pile_up_model': "coates"
}

data_with_pileup = np.copy(data)
IMP.bff.decay_add_pile_up_to_model(
    model=data_with_pileup,
    data=data,
    **pile_up_parameter
)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.semilogy(data_with_pileup, label='with pileup')
ax.semilogy(data, label='no pileup')
ax.set_ylim(10, 100)
ax.legend()

