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
   decays/index
   single-molecule/index
   _imaging/index
   _programming/index

|

.. note:: **Doctest Mode**

   The code-examples in the above tutorials are written in a
   *python-console* format. If you wish to easily execute these examples
   in **IPython**, use::

	%doctest_mode

   in the IPython-console. You can then simply copy and paste the examples
   directly into IPython without having to worry about removing the **>>>**
   manually.