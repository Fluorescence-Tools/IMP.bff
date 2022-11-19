import os
import typing
import pathlib

import RMF
import IMP.atom


def read_angle_file(
        hier: IMP.atom.Hierarchy,
        flex_dict: dict
) -> typing.Tuple[
    typing.List[IMP.atom.Residue],
    typing.List[IMP.atom.Bond]
]:
    flexible_residues = list()
    bonds = list()
    for fr in flex_dict["Flexible residues"]:
        chain_id = fr['chain_identifier']
        residue_id = fr['residue_seq_number']
        sel = IMP.atom.Selection(
            hierarchy=hier,
            chain_ids=[chain_id],
            residue_indexes=[residue_id]
        )
        flexible_residues.append(
            sel.get_selected_particles(False)[0]
        )
    for bnd in flex_dict["Bonds"]:
        a1 = bnd[0]
        a2 = bnd[1]
        sel1 = IMP.atom.Selection(
            hierarchy=hier,
            chain_ids=[a1['chain_identifier']],
            residue_indexes=[a1['residue_seq_number']],
            atom_type=IMP.atom.AtomType(a1['atom_name'])
        )
        sel2 = IMP.atom.Selection(
            hierarchy=hier,
            chain_ids=[a2['chain_identifier']],
            residue_indexes=[a2['residue_seq_number']],
            atom_type=IMP.atom.AtomType(a2['atom_name'])
        )
        p1 = sel1.get_selected_particles()[0]
        p2 = sel2.get_selected_particles()[0]
        b1 = IMP.atom.Bonded(p1)
        b2 = IMP.atom.Bonded(p2)
        bonds.append(
            IMP.atom.create_bond(b1, b2, IMP.atom.Bond.SINGLE)
        )
    return flexible_residues, bonds


class WriteRMFFrame(IMP.OptimizerState):

    def __init__(
            self,
            filename,
            root_hier: IMP.atom.Hierarchy,
            restraints: typing.List[IMP.Restraint],
            name: str = "WriteRMFFrame"
    ):
        model = root_hier.get_model()
        super().__init__(model, name)
        self.restraints = restraints
        fileName, fileExtension = os.path.splitext(filename)
        fn = pathlib.Path(fileName + ".0.rmf3")
        if pathlib.Path(fn).exists():
            for i in range(100):
                fn = pathlib.Path(fileName + ".{:d}.rmf3".format(i))
                if not fn.exists():
                    break
        self._rmf_filename = str(fn)
        self._rh = RMF.create_rmf_file(str(fn))
        IMP.rmf.add_hierarchies(self._rh, [root_hier])
        IMP.rmf.add_restraints(self._rh, restraints)
        IMP.rmf.save_frame(self._rh)

    def do_update(self, arg0):
        #print(*[r.evaluate(False) for r in self.restraints], sep="\t")
        IMP.rmf.save_frame(self._rh)


class WritePDBFrame(IMP.OptimizerState):

    def __init__(
            self,
            filename,
            root_hier: IMP.atom.Hierarchy,
            restraints: typing.List[IMP.Restraint],
            restraint_filename: str = None,
            name: str = "WriteDCDFrame",
            output_objects: typing.List = None,
            multi_state: bool = False
    ):
        model = root_hier.get_model()
        super().__init__(model, name)

        import os
        self._pdb_basename = filename
        if restraint_filename is None:
            restraint_filename = os.path.splitext(filename)[0] + ".rst.txt"
        self._restraint_filename = restraint_filename
        self.restraints = restraints
        self._hier = root_hier
        self.frame = 0
        self.output_objects = output_objects
        self.multi_state = multi_state

    def do_update(self, arg0):
        with open(self._restraint_filename, "a+") as fp:
            fp.write("%s\t" % self.frame)
            fp.write("\t".join(["{:.3f}".format(r.evaluate(False)) for r in self.restraints]))
            fp.write("\t")
            if isinstance(self.output_objects, list):
                for obj in self.output_objects:
                    fp.write(str(obj) + "\t")
            fp.write("\n")
        lead = os.path.splitext(self._pdb_basename)[0]
        if not self.multi_state:
            hiers = [self._hier]
        else:
            hiers = self._hier.get_children()
        for i, hier in enumerate(hiers):
            out_fn = lead + "_state_" + str(i) + "_" + "{:04d}".format(self.frame) + ".pdb"
            IMP.atom.write_pdb(hier, out=out_fn)
        self.frame += 1
