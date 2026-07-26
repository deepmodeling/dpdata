from __future__ import annotations

import dpdata.formats.rdkit.utils
from dpdata.format import Format


@Format.register("mol")
@Format.register("mol_file")
class MolFormat(Format):
    """MDL Molfile containing one molecular graph and its conformers.

    `RDKit <https://www.rdkit.org/>`_ is a collection of cheminformatics and
    machine-learning tools.

    Reading and writing requires RDKit. Bond orders and formal charges are
    preserved through :class:`dpdata.BondOrderSystem` rather than the regular
    System classes.
    """

    def from_bond_order_system(self, file_name, **kwargs):
        """Load an MDL Molfile as an RDKit molecule.

        Parameters
        ----------
        file_name : str or os.PathLike
            Input ``.mol`` file.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        rdkit.Chem.Mol
            Molecule with explicit hydrogens retained and sanitization deferred
            to :class:`dpdata.BondOrderSystem`.
        """
        import rdkit.Chem

        return rdkit.Chem.MolFromMolFile(file_name, sanitize=False, removeHs=False)

    def to_bond_order_system(self, data, mol, file_name, frame_idx=0, **kwargs):
        """Write one conformer to an MDL Molfile.

        Parameters
        ----------
        data : dict
            BondOrderSystem data.
        mol : rdkit.Chem.Mol
            RDKit molecule carrying the bond graph and conformers.
        file_name : str or os.PathLike
            Destination ``.mol`` file.
        frame_idx : int, default=0
            Conformer/frame index to write.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.
        """
        import rdkit.Chem

        assert frame_idx < mol.GetNumConformers()
        rdkit.Chem.MolToMolFile(mol, file_name, confId=frame_idx)


@Format.register("sdf")
@Format.register("sdf_file")
class SdfFormat(Format):
    """Structure-data file (SDF) containing one or more conformers.

    `RDKit <https://www.rdkit.org/>`_ provides cheminformatics capabilities
    for SDF reading and writing.

    All records must describe the same molecular topology so they can be
    represented as conformers of one :class:`dpdata.BondOrderSystem`. Reading
    and writing requires RDKit.
    """

    def from_bond_order_system(self, file_name, **kwargs):
        """Load same-topology SDF records as one RDKit molecule.

        Parameters
        ----------
        file_name : str or os.PathLike
            Input ``.sdf`` file. All records must share a topology.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        rdkit.Chem.Mol
            Molecule whose conformers correspond to the SDF records.
        """
        import rdkit.Chem

        mols = [
            m
            for m in rdkit.Chem.SDMolSupplier(file_name, sanitize=False, removeHs=False)
        ]
        if len(mols) > 1:
            mol = dpdata.formats.rdkit.utils.combine_molecules(mols)
        else:
            mol = mols[0]
        return mol

    def to_bond_order_system(self, data, mol, file_name, frame_idx=-1, **kwargs):
        """Write conformers to an SDF file.

        Parameters
        ----------
        data : dict
            BondOrderSystem data.
        mol : rdkit.Chem.Mol
            RDKit molecule carrying the bond graph and conformers.
        file_name : str or os.PathLike
            Destination ``.sdf`` file.
        frame_idx : int, default=-1
            Conformer to write. ``-1`` writes every conformer as a separate
            SDF record.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.
        """
        import rdkit.Chem

        sdf_writer = rdkit.Chem.SDWriter(file_name)
        if frame_idx == -1:
            for ii in range(mol.GetNumConformers()):
                sdf_writer.write(mol, confId=ii)
        else:
            assert frame_idx < mol.GetNumConformers()
            sdf_writer.write(mol, confId=frame_idx)
        sdf_writer.close()
