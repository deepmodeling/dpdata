import numpy as np

from dpdata.periodic_table import ELEMENTS
from dpdata.unit import EnergyConversion

kcal2ev = EnergyConversion("kcal_mol", "eV").value()

START = 0
READ_CHARGE = 2
READ_COORDS_START = 3
READ_COORDS = 6
READ_FORCES = 7


def parse_sqm_out(fname):
    """Read atom symbols, charges and coordinates from ambertools sqm.out file."""
    atom_symbols = []
    coords = []
    charges = []
    forces = []
    energies = []

    with open(fname) as f:
        flag = START
        for line in f:
            if line.startswith(" Total SCF energy"):
                energy = float(line.strip().split()[-2])
                energies = [energy]
            elif line.startswith("  Atom    Element       Mulliken Charge"):
                flag = READ_CHARGE
                charges = []
            elif line.startswith(" Total Mulliken Charge"):
                flag = START
            elif line.startswith(" Final Structure"):
                flag = READ_COORDS_START
                coords = []
            elif line.startswith("QMMM: Forces on QM atoms"):
                flag = READ_FORCES
                forces = []
            elif flag == READ_CHARGE:
                ls = line.strip().split()
                atom_symbols.append(ls[-2])
                charges.append(float(ls[-1]))
            elif READ_COORDS_START <= flag < READ_COORDS:
                flag += 1
            elif flag == READ_COORDS:
                coords.append([float(x) for x in line.strip().split()[-3:]])
                if len(coords) == len(charges):
                    flag = START
            elif flag == READ_FORCES:
                ll = line.strip()
                if not ll.startswith("QMMM: Atm "):
                    flag = START
                    continue
                forces.append([float(ll[-60:-40]), float(ll[-40:-20]), float(ll[-20:])])
                if len(forces) == len(charges):
                    flag = START

    data = {}
    atom_names, data["atom_types"], atom_numbs = np.unique(
        atom_symbols, return_inverse=True, return_counts=True
    )
    data["charges"] = np.array(charges)
    data["atom_names"] = list(atom_names)
    data["atom_numbs"] = list(atom_numbs)
    data["orig"] = np.array([0, 0, 0])
    data["cells"] = np.array(
        [[[100.0, 0.0, 0.0], [0.0, 100.0, 0.0], [0.0, 0.0, 100.0]]]
    )
    data["nopbc"] = True
    data["coords"] = np.array([coords])

    energies = np.array(energies)
    forces = -np.array([forces], dtype=np.float64) * kcal2ev
    if len(forces) > 0:
        data["energies"] = energies
        data["forces"] = forces

    return data


def make_sqm_in(data, fname=None, frame_idx=0, **kwargs):
    symbols = [data["atom_names"][ii] for ii in data["atom_types"]]
    atomic_numbers = [ELEMENTS.index(ss) + 1 for ss in symbols]
    charge = kwargs.get("charge", 0)

    # multiplicity
    mult = kwargs.get("mult", 1)
    if mult != 1:
        raise RuntimeError("Multiplicity is not 1, which is not supported by sqm")

    maxcyc = kwargs.get("maxcyc", 0)  # 0 represents a single-point calculation
    theory = kwargs.get("qm_theory", "DFTB3")
    ret = "Run semi-emperical minimization\n"
    ret += " &qmmm\n"
    ret += f"     qm_theory='{theory}'\n"
    ret += f"     qmcharge={charge}\n"
    ret += f"     maxcyc={maxcyc}\n"
    ret += "     verbosity=4\n"
    ret += " /\n"
    for ii in range(len(data["atom_types"])):
        ret += "{:>4s}{:>6s}{:>16s}{:>16s}{:>16s}\n".format(
            str(atomic_numbers[ii]),
            str(symbols[ii]),
            f"{data['coords'][frame_idx][ii, 0]:.6f}",
            f"{data['coords'][frame_idx][ii, 1]:.6f}",
            f"{data['coords'][frame_idx][ii, 2]:.6f}",
        )
    if fname is not None:
        with open(fname, "w") as fp:
            fp.write(ret)
    return ret
