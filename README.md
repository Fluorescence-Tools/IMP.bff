\brief Bayesian Fluorescence Framework

# imp_bff {#imp_bff}

Bayesian Fluorescence Framework (BFF) processes and analyzes 
fluorescence data. BFF contains functions to computed label 
distributions and forward models. Among other data the program 
imp_bff samples can sample biomolecular conformations restrained 
by experimental data by  simulating fluorophore distributions 
around attachment sites and  by comparing simulated observables 
to experimental data.


## Inter-label distance score usage:

First, import the module:

```
import IMP.bff
import IMP.bff.restraints
```

Then, select the "score set", i.e., a set of distances that are used for 
score calculation from a FPS.JSON file:

```
fps_json_fn = str(root_dir / "screening.fps.json")
score_set = "inter"
```

Finally, create the restraint and add it to the model.

```
fret_restraint = IMP.bff.restraints.AVNetworkRestraintWrapper(
    hier, fps_json_fn,
    mean_position_restraint=True,
    score_set=score_set
)
fret_restraint.add_to_model()
output_objects.append(fret_restraint)
```


# Info

_Author(s)_: Thomas-Otavio Peulen

_Maintainer_: `tpeulen`

_License_: [LGPL](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html)
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Publications_:
- None
