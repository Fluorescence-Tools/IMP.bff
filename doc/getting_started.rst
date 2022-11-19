.. This helps define the TOC ordering for "about us" sections. Particularly
   useful for PDF output as this section is not linked from elsewhere.

.. Places global toc into the sidebar

:globalsidebartoc: True

.. _preface_menu:

.. include:: includes/big_toc_css.rst
.. include:: tune_toc.rst

.. _getting_started:
Getting Started
===============
``IMP.bff`` is an open source framework for Bayesian analysis of fluorescence
data. BFF provides various tools for processing and data preprocessing.

IMP.bff is not an application in the traditional sense: there is no program called “IMP.bff”
that you run. IMP.bff is a collection of libraries and different programs written in C++ and
the Python programming language. Those libraries and programs can easily be chained together to
create tools and pipelines for integrative modeling with fluorescence information with the
integrative modeling platform (IMP). Don’t worry, you don’t need to know anything about Python
programming to use existing simulation pipelines. You can use the existing application examples
as templates.

On the other hand, if you don’t mind doing a little programming, this approach gives you enormous
power and flexibility. Your script has complete access to the entire IMP and IMP.bff application
programming interface (API), as well as the full power of the Python language and libraries.
You have complete control over every detail of the integrative modeling process, from defining
the representation of your system to analyzing the results.

The purpose of this guide is to illustrate some of the main ``IMP.bff`` features. It assumes a very
basic working knowledge on the integrative modeling platform (IMP), fluorescence spectroscopy
(decay analysis, correlation spectroscopy, etc.). Please refer to our :ref:`installation instructions
<installation-instructions>` for installing ``IMP.bff``.


Accessible volume simulations
-----------------------------
``IMP.bff`` implements a IMP decorator for accessible volumes (AVs). An AV is the
sterically allowed conformational space of a label around its attachment site
(source).

 .. code-block:: python

    import IMP
    import IMP.atom
    import IMP.bff

    m = IMP.Model()
    source = IMP.Particle(m)
    IMP.core.XYZ.setup_particle(source, (0, 0, 2))

    av_p = IMP.Particle(m)
    av_parameter = {
        "linker_length": 20.0,
        "radii": (3.5, 0.0, 0.0),
        "linker_width": 0.5,
        "allowed_sphere_radius": 2.0,
        "contact_volume_thickness": 0.0,
        "contact_volume_trapped_fraction": -1,
        "simulation_grid_resolution": 0.5
    }
    IMP.bff.AV.do_setup_particle(m, av_p, source, **av_parameter)
    av1 = IMP.bff.AV(m, av_p)

Distance distributions between AVs can be computed by sampling over two AVs

 .. code-block:: python

    import numpy as np
    rda_start, rda_stop, n_bins = 0, 100, 128
    n_samples = 10000
    rda = np.linspace(rda_start, rda_stop, n_bins)
    p_rda = IMP.bff.av_distance_distribution(av1, av2, rda, n_samples=n_samples)

Such distance distributions can be used to compute experimental observables.


Restrained integrative modeling
-------------------------------
A set of simulated labeled can be used as restraint for integrative modeling.
To incoporate a restain for integratove modeling, first, import the module:

 .. code-block:: python

    import IMP.bff
    import IMP.bff.restraints

Then, define a fps.json file that contain experimental information:

 .. code-block:: json

    {
      "Distances": {
        "5-44_C3": {
          "error_neg": 3.6,
          "error_pos": 3.6,
          "distance": 25.6,
          "position1_name": "5D",
          "position2_name": "44A",
          "Forster_radius": 52,
          "distance_type": "RDAMean"
        },
        "5-132_C3": {
          "error_neg": 5.9,
          "error_pos": 5.9,
          "distance": 35.8,
          "position1_name": "5D",
          "position2_name": "132A",
          "Forster_radius": 52,
          "distance_type": "RDAMean"
        }
      },
      "Positions": {
        "5D": {
            "allowed_sphere_radius": 5,
            "atom_name": "CB",
            "linker_length": 20,
            "linker_width": 2.5,
            "radius1": 3.5,
            "residue_seq_number": 5,
            "simulation_grid_resolution": 2.0,
            "simulation_type": "AV1"
        },
        "132A": {
            "atom_name": "CB",
            "contact_volume_thickness": 3,
            "linker_length": 22,
            "linker_width": 3.5,
            "radius1": 3.5,
            "residue_seq_number": 132,
            "simulation_grid_resolution": 2.0,
            "simulation_type": "AV1"
        },
        "44A": {
            "atom_name": "CB",
            "contact_volume_thickness": 3,
            "linker_length": 22,
            "linker_width": 4.5,
            "radius1": 3.5,
            "residue_seq_number": 44,
            "simulation_grid_resolution": 2.0,
            "simulation_type": "AV1"
        }
      },
      "χ²": {
        "chi2_C3_all": {
          "distances": [
            "5-132_C3",
            "5-44_C3"
          ]
        },
        "chi2_C3": {
          "distances": [
            "5-132_C2"
          ]
      }
    }

In your modeling script import the module:

 .. code-block:: python

    import IMP.bff
    import IMP.bff.restraints


Finally, select the "score set", i.e., a set of distances that are used for
score calculation, create the restraint and add it to the model.

 .. code-block:: python

    fps_json_fn = str(root_dir / "screening.fps.json")
    score_set = "chi2_C3_all"
    fret_restraint = IMP.bff.restraints.AVNetworkRestraintWrapper(
        hier, fps_json_fn,
        mean_position_restraint=True,
        score_set=score_set
    )
    fret_restraint.add_to_model()


The restraint has two modes of operation. Either the AVs are computed once
and the label mean positions are attached to rigid bodies and the distance
between the mean positions are used for scoring (mean_position_restraint) or
the AVs are recomputed at every iteration. The later is computationally more
expensive (but more accurate). More features of ``IMP.bff`` are described the
examples.

Next steps
----------
We have briefly covered how ``IMP.bff`` can be used for restrained simulations.
Please refer to our :ref:`user_guide` for details on all the tools that we provide.
You can also find an exhaustive list of the public API in the
:ref:`api_ref` or the numerous :ref:`examples` that
illustrate the use of ``IMP.bff`` in many different contexts.

The :ref:`tutorials <tutorial_menu>` also contain additional learning
resources.
