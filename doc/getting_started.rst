Getting Started
===============
The purpose of this guide is to illustrate some of the main features that
``IMP.bff`` provides. It assumes a very basic working knowledge on the
integrative modeling platform (IMP), fluorescence spectroscopy (decay analysis,
correlation spectroscopy, etc.). Please refer to our :ref:`installation instructions
<installation-instructions>` for installing ``IMP.bff``.

``IMP.bff`` is an open source framework for Bayesian analysis of fluorescence
data. BFF provides various tools for processing and data preprocessing.

AV decorator
------------

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

AVs are computed using path search algorithms (see examples llllll).

PMI restraints
--------------


First, import the module:

 .. code-block:: python

    import IMP.bff
    import IMP.bff.restraints

Then, select the "score set", i.e., a set of distances that are used for
score calculation from a FPS.JSON file:

 .. code-block:: python

    fps_json_fn = str(root_dir / "screening.fps.json")
    score_set = "inter"

Finally, create the restraint and add it to the model.

 .. code-block:: python

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
expensive (but more accurate).

Next steps
----------
We have briefly covered how ``IMP.bff`` can be used for restrained simuations.
Please refer to our :ref:`user_guide` for details on all the tools that we provide.
You can also find an exhaustive list of the public API in the
:ref:`api_ref` or the numerous :ref:`examples` that
illustrate the use of ``IMP.bff`` in many different contexts.

The :ref:`tutorials <tutorial_menu>` also contain additional learning
resources.
