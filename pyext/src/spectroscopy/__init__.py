"""Transitional package -- its contents move out in PRD-113 stage 4.

``decay.py`` is gone: it wrapped the C++ ``Decay*`` classes to convolve a decay
with an IRF, apply pileup and score the result against measured counts. That is
the experiment, and ``IMP.bff`` emits experiment-neutral quantities -- lifetime
spectra and rate constants, unconvolved. The whole C++ instrument layer went
with it (see the module log).

``kappa2.py`` stays for now and is the opposite kind of thing: a forward model
producing orientation-factor distributions, wobbling-in-a-cone averages and
order parameters from structure. It moves to ``IMP.bff.photophysics`` in stage 4,
joining ``fret/kappa2.py`` so that "kappa squared" has one home rather than two.
"""
