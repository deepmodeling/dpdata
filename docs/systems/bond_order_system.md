
## BondOrderSystem
A new class {class}`BondOrderSystem` which inherits from class {class}`System` is introduced in dpdata. This new class contains information of chemical bonds and formal charges (stored in `BondOrderSystem.data['bonds']`, `BondOrderSystem.data['formal_charges']`). Now BondOrderSystem can only read from .mol/.sdf formats, because of its dependency on rdkit (which means rdkit must be installed if you want to use this function). Other formats, such as pdb, must be converted to .mol/.sdf format (maybe with software like open babel).
```python
import dpdata

system_1 = dpdata.BondOrderSystem(
    "tests/bond_order/CH3OH.mol", fmt="mol"
)  # read from .mol file
system_2 = dpdata.BondOrderSystem(
    "tests/bond_order/methane.sdf", fmt="sdf"
)  # read from .sdf file
```
In sdf file, all molecules must be of the same topology (i.e. conformers of the same molecular configuration).
`BondOrderSystem` also supports initialize from a {class}`rdkit.Chem.rdchem.Mol` object directly.
```python
from rdkit import Chem
from rdkit.Chem import AllChem
import dpdata

mol = Chem.MolFromSmiles("CC")
mol = Chem.AddHs(mol)
AllChem.EmbedMultipleConfs(mol, 10)
system = dpdata.BondOrderSystem(rdkit_mol=mol)
```

### Bond Order Assignment
The {class}`BondOrderSystem` implements a more robust sanitize procedure for rdkit Mol, as defined in {class}`dpdata.rdkit.santizie.Sanitizer`. This class defines 3 level of sanitization process by: low, medium and high. (default is medium).
+ low: use `rdkit.Chem.SanitizeMol()` function to sanitize molecule.
+ medium: before using rdkit, the programm will first assign formal charge of each atom to avoid inappropriate valence exceptions. However, this mode requires the rightness of the bond order information in the given molecule.
+ high: the program will try to fix inappropriate bond orders in aromatic hetreocycles, phosphate, sulfate, carboxyl, nitro, nitrine, guanidine groups. If this procedure fails to sanitize the given molecule, the program will then try to call `obabel` to pre-process the mol and repeat the sanitization procedure. **That is to say, if you wan't to use this level of sanitization, please ensure `obabel` is installed in the environment.**
According to our test, our sanitization procedure can successfully read 4852 small molecules in the PDBBind-refined-set. It is necessary to point out that the in the molecule file (mol/sdf), the number of explicit hydrogens has to be correct. Thus, we recommend to use
 `obabel xxx -O xxx -h` to pre-process the file. The reason why we do not implement this hydrogen-adding procedure in dpdata is that we can not ensure its correctness.

```python
import dpdata

for sdf_file in glob.glob("bond_order/refined-set-ligands/obabel/*sdf"):
    syst = dpdata.BondOrderSystem(sdf_file, sanitize_level="high", verbose=False)
```
### Formal Charge Assignment
BondOrderSystem implement a method to assign formal charge for each atom based on the 8-electron rule (see below). Note that it only supports common elements in bio-system: B,C,N,O,P,S,As
```python
import dpdata

syst = dpdata.BondOrderSystem("tests/bond_order/CH3NH3+.mol", fmt="mol")
print(syst.get_formal_charges())  # return the formal charge on each atom
print(syst.get_charge())  # return the total charge of the system
```

If a valence of 3 is detected on carbon, the formal charge will be assigned to -1. Because for most cases (in alkynyl anion, isonitrile, cyclopentadienyl anion), the formal charge on 3-valence carbon is -1, and this is also consisent with the 8-electron rule.
