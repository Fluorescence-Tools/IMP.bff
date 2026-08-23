"""PRD-113 stage 2: a labelling site and a quencher are things, not dicts.

Before this, a site was an untyped ``source_info`` dict threaded through the AV
builder, the quenching model and every benchmark -- and that dict mixed *where
the dye is attached* with *how its accessible volume is computed*. Quenching
parameters were worse: ``kQ`` and ``rC`` were looked up at a call site and
passed as bare floats through six layers, which is why PRD-111 had to build a
``quencher_table(kQ_scale, rC)`` helper inside a benchmark to vary them at all.

Names follow the FLR dictionaries where an item exists (checked against
``../mmfdb/src/mmfdb/data/*.dic``); quenching has no item anywhere in the stack.
"""

import pytest

from IMP.bff import find_dye
from IMP.bff import (
    Label, Quencher, PETParameters, pet_quenching_reference,
    reference_quenchers, reference_pet_parameters, REFERENCE_DYE,
    FLUOROPHORE_TYPES,
)
from IMP.bff import LABEL_FLRCIF_ITEMS as SITE_ITEMS
from IMP.bff import QUENCHER_FLRCIF_ITEMS as QUENCHER_ITEMS


class TestLabel:

    def test_a_label_is_a_position_plus_a_dye(self):
        label = Label(asym_id="A", seq_id=132, atom_id="CB",
                      dye=find_dye("AlexaFluor 488"), fluorophore_type="donor")
        assert label.key == ("A", 132, "CB")
        assert label.dye.name == "AlexaFluor488"

    def test_the_fluorophore_type_enum_is_the_dictionary_s(self):
        """``_flr_sample_probe_details.fluorophore_type``, verbatim."""
        assert FLUOROPHORE_TYPES == ("donor", "acceptor", "unspecified")
        with pytest.raises(ValueError, match="fluorophore_type"):
            Label(asym_id="A", seq_id=1, fluorophore_type="quencher")

    def test_it_round_trips_the_fps_dialect(self):
        label = Label(asym_id="E", seq_id=96, atom_id="CB")
        assert Label.from_source_info(label.to_source_info()) == label

    def test_it_carries_only_the_position_half(self):
        """The AV parameters in a source_info dict are *not* a site property.

        ``linker_length`` and ``allowed_sphere_radius`` describe how an
        accessible volume is computed; a rotamer library has neither, and a
        label does not change when the representation does.
        """
        source_info = {
            "chain_identifier": "A", "residue_seq_number": 132, "atom_name": "CB",
            "linker_length": 20.0, "linker_width": 0.5, "radius1": 3.5,
            "allowed_sphere_radius": 2.1, "simulation_grid_resolution": 1.5,
        }
        label = Label.from_source_info(source_info)
        carried = set(label.to_source_info())
        assert carried == {"chain_identifier", "residue_seq_number", "atom_name"}
        for representation_parameter in (
                "linker_length", "radius1", "allowed_sphere_radius",
                "simulation_grid_resolution"):
            assert not hasattr(label, representation_parameter)

    def test_the_position_names_are_the_dictionary_s(self):
        for field, item in SITE_ITEMS.items():
            if item is not None:
                assert item.startswith("_flr_"), f"{field} -> {item}"


class TestQuencher:

    def test_the_redox_atoms_are_not_CB(self):
        """Electron transfer happens at the ring, the thioether or the thiol.

        Stamping the rate on CB puts it up to 4 A from where the chemistry is.
        """
        quenchers = reference_quenchers()
        assert "CB" not in quenchers["TRP"].atom_ids
        assert "NE1" in quenchers["TRP"].atom_ids     # indole nitrogen
        assert quenchers["MET"].atom_ids == ("SD",)   # thioether sulfur
        assert "OH" in quenchers["TYR"].atom_ids      # phenol oxygen

    def test_a_quencher_carries_no_rate(self):
        """PET takes two partners, so kQ cannot belong to the tryptophan.

        A rhodamine, an oxazine and a cyanine see the same tryptophan
        differently -- the rate depends on both redox potentials.
        """
        trp = reference_quenchers()["TRP"]
        assert not hasattr(trp, "rate_constant")
        assert not hasattr(trp, "attenuation_length")

    def test_it_reads_the_one_table_rather_than_restating_it(self):
        reference = pet_quenching_reference()
        params = reference_pet_parameters()
        assert set(params) == set(reference)
        for comp_id, entry in reference.items():
            assert params[comp_id].rate_constant == pytest.approx(entry.kQ)
            assert params[comp_id].contact_distance == pytest.approx(
                entry.contact_distance)

    def test_the_rate_is_the_knob_a_calibration_turns(self):
        """The published kQ are starting values, not constants."""
        base = reference_pet_parameters()["TRP"].rate_constant
        assert reference_pet_parameters(rate_scale=0.5)["TRP"].rate_constant == pytest.approx(0.5 * base)
        assert base * 0.25 == pytest.approx(
            reference_pet_parameters()["TRP"].scaled(0.25).rate_constant)

    def test_transfer_to_another_dye_is_visible(self):
        """The bundled table is for a xanthene; applying it elsewhere is an
        assumption, and it should be possible to see that it was made."""
        native = reference_pet_parameters(REFERENCE_DYE)["TRP"]
        assert not native.is_transferred
        other = reference_pet_parameters("Atto655")["TRP"]
        assert other.is_transferred
        assert other.measured_for == REFERENCE_DYE
        assert other.rate_constant == native.rate_constant

    def test_pet_is_keyed_by_the_pair(self):
        p = PETParameters(dye="Atto655", comp_id="TRP", rate_constant=1.0,
                          contact_distance=5.0)
        assert (p.dye, p.comp_id) == ("Atto655", "TRP")

    def test_a_type_can_be_located_in_a_structure(self):
        trp = reference_quenchers()["TRP"]
        assert trp.is_typed
        here = trp.at("A", 126)
        assert not here.is_typed
        assert (here.asym_id, here.seq_id) == ("A", 126)
        assert here.atom_ids == trp.atom_ids

    def test_quenching_has_no_dictionary_item(self):
        """Recorded so nobody later assumes it was an oversight.

        Checked across all ten .dic files in ../mmfdb/src/mmfdb/data: nothing
        matches "quench". The identifiers do have items; the photophysics does not.
        """
        assert QUENCHER_ITEMS["rate_constant"] is None
        assert QUENCHER_ITEMS["attenuation_length"] is None
        assert QUENCHER_ITEMS["comp_id"] == "_flr_poly_probe_position.comp_id"

    def test_a_negative_rate_is_refused(self):
        with pytest.raises(ValueError, match="rate_constant"):
            PETParameters(dye="d", comp_id="TRP", rate_constant=-1.0,
                          contact_distance=5.0)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
