
 .. Spectroscopy:

************
Spectroscopy
************


Forster resonance energy transfer
################################


Orientation factor
##################


Fluorescence decays
###################

Phasor
******
Changing the data representation from the classical time delay
histogram to the phasor representation provides a global view
of the fluorescence decay at each pixel of an image. In the phasor
representation reveal the presence of different in a a fluorescence
decay. The analysis of time-resolve fluorescence data in the phasor space is done
by observing clustering of populations in the phasor plot rather
than by fitting the fluorescence decay using exponentials. The
analysis is instantaneous since is not based on calculations or nonlinear
fitting. Thus, the phasor approach can simplify the way data
are analyzed and makes FLIM technique accessible to the nonexpert
in spectroscopy and data analysis :cite:`DIGMAN2008L14`.

The phasor transformation can use data collected in either the time or the
frequency domain. In the time-correlated single-photon counting (TCSPC)
technique, the histogram of the photon arrival time at each pixel of the image
is transformed to Fourier space, and the phasor coordinates are calculated
using the following relations:

.. math::

    g_{i,j}(\omega) = \int_{0}^{T} I(t) \cdot \cos(n\omega t) dt / \int_{0}^{T} I(t) dt

    s_{i,j}(\omega) = \int_{0}^{T} I(t) \cdot \sin(n\omega t) dt / \int_{0}^{T} I(t) dt

in which gi,j(ω) and si,j(ω) are the x and y coordinates of the phasor plot,
n and ω are the harmonic frequency and the angular frequency of excitation,
respectively, and T is the repeat frequency of the acquisition :cite:`ranjit2018fit`.

