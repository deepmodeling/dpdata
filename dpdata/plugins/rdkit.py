from __future__ import annotations

import dpdata.rdkit.utils
from dpdata.format import Format


@Format.register("mol")
@Format.register("mol_file")
class MolFormat(Format):
    def from_bond_order_system(self, file_name, **kwargs):
        import rdkit.Chem

        return rdkit.Chem.MolFromMolFile(file_name, sanitize=False, removeHs=False)

    def to_bond_order_system(self, data, mol, file_name, frame_idx=0, **kwargs):
        import rdkit.Chem

        assert frame_idx < mol.GetNumConformers()
        rdkit.Chem.MolToMolFile(mol, file_name, confId=frame_idx)


@Format.register("sdf")
@Format.register("sdf_file")
class SdfFormat(Format):
    def from_bond_order_system(self, file_name, **kwargs):
        """Note that it requires all molecules in .sdf file must be of the same topology."""
        import rdkit.Chem

        mols = [
            m
            for m in rdkit.Chem.SDMolSupplier(file_name, sanitize=False, removeHs=False)
        ]
        if len(mols) > 1:
            mol = dpdata.rdkit.utils.combine_molecules(mols)
        else:
            mol = mols[0]
        return mol

    def to_bond_order_system(self, data, mol, file_name, frame_idx=-1, **kwargs):
        import rdkit.Chem

        sdf_writer = rdkit.Chem.SDWriter(file_name)
        if frame_idx == -1:
            for ii in range(mol.GetNumConformers()):
                sdf_writer.write(mol, confId=ii)
        else:
            assert frame_idx < mol.GetNumConformers()
            sdf_writer.write(mol, confId=frame_idx)
        sdf_writer.close()
