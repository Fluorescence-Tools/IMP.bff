.. Places global toc into the sidebar

:globalsidebartoc: True

.. _tutorial_menu:


.. include:: ../includes/big_toc_css.rst
.. include:: ../tune_toc.rst

=================
User manual
=================

|

.. toctree::
   :maxdepth: 2
   :glob:

   _concepts/index
   structure/index
   _intensity/index
   fcs/index
   single-molecule/index
   _imaging/index

|

.. note:: **Time-resolved decays are not here.**

   ``IMP.bff`` emits experiment-neutral quantities -- lifetime spectra, rate
   constants, :math:`\kappa^2` distributions, distances. Everything that turns
   one into a measured histogram -- convolution with an IRF, pile-up,
   linearisation, counting statistics and decay fitting -- lives in
   ``tttrlib`` (``modules/spectroscopy/decay``: the ``fconv`` family,
   ``BlindIRF``, ``DecayFit23``--``26``, ``MaxEntTcspc``, ``DecayStatistics``).
   A forward model built here is handed to ``tttrlib`` to be compared with
   photons.

|

.. note:: **Doctest Mode**

   The code-examples in the above tutorials are written in a
   *python-console* format. If you wish to easily execute these examples
   in **IPython**, use::

	%doctest_mode

   in the IPython-console. You can then simply copy and paste the examples
   directly into IPython without having to worry about removing the **>>>**
   manually.