
 .. BiomolecularStructure:

===================
Structural modeling
===================

Accessible volume restraint
---------------------------


```python
import IMP.bff
import IMP.bff.restraints
```

Then, select the "score set", i.e., a set of distances that are used for
score calculation from a FPS.JSON file:

```python
fps_json_fn = str(root_dir / "screening.fps.json")
score_set = "inter"
```

Finally, create the restraint and add it to the model.

```python
fret_restraint = IMP.bff.restraints.AVNetworkRestraintWrapper(
    hier, fps_json_fn,
    mean_position_restraint=True,
    score_set=score_set
)
fret_restraint.add_to_model()
output_objects.append(fret_restraint)
```

