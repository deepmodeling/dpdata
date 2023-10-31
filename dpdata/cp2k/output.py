# %%
import math
import re
from collections import OrderedDict

import numpy as np

from ..unit import (
    EnergyConversion,
    ForceConversion,
    LengthConversion,
    PressureConversion,
)
from .cell import cell_to_low_triangle

AU_TO_ANG = LengthConversion("bohr", "angstrom").value()
AU_TO_EV = EnergyConversion("hartree", "eV").value()
AU_TO_EV_EVERY_ANG = ForceConversion("hartree/bohr", "eV/angstrom").value()
delimiter_patterns = []
delimiter_p1 = re.compile(r"^ \* GO CP2K GO! \*+")
delimiter_p2 = re.compile(r"^ \*+")
delimiter_patterns.append(delimiter_p1)
delimiter_patterns.append(delimiter_p2)
avail_patterns = []
avail_patterns.append(re.compile(r"^ INITIAL POTENTIAL ENERGY"))
avail_patterns.append(re.compile(r"^ ENSEMBLE TYPE"))


class Cp2kSystems:
    """deal with cp2k outputfile."""

    def __init__(self, log_file_name, xyz_file_name, restart=False):
        self.log_file_object = open(log_file_name)
        self.xyz_file_object = open(xyz_file_name)
        self.log_block_generator = self.get_log_block_generator()
        self.xyz_block_generator = self.get_xyz_block_generator()
        self.restart_flag = restart

        self.cell = None
        self.print_level = None

        self.atomic_kinds = None

        if self.restart_flag:
            self.handle_single_log_frame(next(self.log_block_generator))

    def __del__(self):
        self.log_file_object.close()
        self.xyz_file_object.close()

    def __iter__(self):
        return self

    def __next__(self):
        info_dict = {}
        log_info_dict = self.handle_single_log_frame(next(self.log_block_generator))
        # print(log_info_dict)
        xyz_info_dict = self.handle_single_xyz_frame(next(self.xyz_block_generator))
        # eq1 = [v1==v2 for v1,v2 in zip(log_info_dict['atom_numbs'], xyz_info_dict['atom_numbs'])]
        # eq2 = [v1==v2 for v1,v2 in zip(log_info_dict['atom_names'], xyz_info_dict['atom_names'])]
        # eq3 = [v1==v2 for v1,v2 in zip(log_info_dict['atom_types'], xyz_info_dict['atom_types'])]
        # assert all(eq1), (log_info_dict,xyz_info_dict,'There may be errors in the file. If it is a restart task; use restart=True')
        # assert all(eq2), (log_info_dict,xyz_info_dict,'There may be errors in the file. If it is a restart task; use restart=True')
        # assert all(eq3), (log_info_dict,xyz_info_dict,'There may be errors in the file. If it is a restart task; use restart=True')
        assert math.isclose(
            log_info_dict["energies"], xyz_info_dict["energies"], abs_tol=1.0e-6
        ), (
            log_info_dict["energies"],
            xyz_info_dict["energies"],
            "There may be errors in the file",
        )
        info_dict.update(log_info_dict)
        info_dict.update(xyz_info_dict)
        return info_dict

    def get_log_block_generator(self):
        lines = []
        delimiter_flag = False
        yield_flag = False
        while True:
            line = self.log_file_object.readline()
            if line:
                lines.append(line)
                if any(p.match(line) for p in delimiter_patterns):
                    if delimiter_flag is True:
                        yield_flag = True
                        yield lines
                        lines = []
                        delimiter_flag = False
                    else:
                        line = self.log_file_object.readline()
                        lines.append(line)
                        if any(p.match(line) for p in avail_patterns):
                            delimiter_flag = True
            else:
                if not yield_flag:
                    raise StopIteration("None of the delimiter patterns are matched")
                break
        if delimiter_flag is True:
            raise RuntimeError("This file lacks some content, please check")

    def get_xyz_block_generator(self):
        p3 = re.compile(r"^\s*(\d+)\s*")
        yield_flag = False
        while True:
            line = self.xyz_file_object.readline()
            if not line:
                if not yield_flag:
                    raise StopIteration("None of the xyz patterns are matched")
                break
            if p3.match(line):
                yield_flag = True
                atom_num = int(p3.match(line).group(1))
                lines = []
                lines.append(line)
                for ii in range(atom_num + 1):
                    lines.append(self.xyz_file_object.readline())
                if not lines[-1]:
                    raise RuntimeError(
                        "this xyz file may lack of lines, should be {};lines:{}".format(
                            atom_num + 2, lines
                        )
                    )
                yield lines

    def handle_single_log_frame(self, lines):
        info_dict = {}
        energy_pattern_1 = re.compile(
            r" INITIAL POTENTIAL ENERGY\[hartree\]\s+=\s+(?P<number>\S+)"
        )
        #  CONSERVED QUANTITY [hartree] =                              -0.279168013085E+04
        energy_pattern_2 = re.compile(
            r" POTENTIAL ENERGY\[hartree\]\s+=\s+(?P<number>\S+)"
        )
        energy = None
        cell_length_pattern = re.compile(
            r" (INITIAL ){0,1}CELL LNTHS\[bohr\]\s+=\s+(?P<A>\S+)\s+(?P<B>\S+)\s+(?P<C>\S+)"
        )
        cell_angle_pattern = re.compile(
            r" (INITIAL ){0,1}CELL ANGLS\[deg\]\s+=\s+(?P<alpha>\S+)\s+(?P<beta>\S+)\s+(?P<gamma>\S+)"
        )
        cell_A, cell_B, cell_C = (
            0,
            0,
            0,
        )
        cell_alpha, cell_beta, cell_gamma = (
            0,
            0,
            0,
        )
        cell_a_pattern = re.compile(
            r" CELL\| Vector a \[angstrom\]:\s+(?P<ax>\S+)\s+(?P<ay>\S+)\s+(?P<az>\S+)"
        )
        cell_b_pattern = re.compile(
            r" CELL\| Vector b \[angstrom\]:\s+(?P<bx>\S+)\s+(?P<by>\S+)\s+(?P<bz>\S+)"
        )
        cell_c_pattern = re.compile(
            r" CELL\| Vector c \[angstrom\]:\s+(?P<cx>\S+)\s+(?P<cy>\S+)\s+(?P<cz>\S+)"
        )
        force_start_pattern = re.compile(r" ATOMIC FORCES in")
        force_flag = False
        force_end_pattern = re.compile(r" SUM OF ATOMIC FORCES")
        force_lines = []
        cell_flag = 0
        print_level_pattern = re.compile(
            r" GLOBAL\| Global print level\s+(?P<print_level>\S+)"
        )
        print_level_flag = 0
        atomic_kinds_pattern = re.compile(r"\s+\d+\. Atomic kind:\s+(?P<akind>\S+)")
        atomic_kinds = []
        stress_sign = "STRESS"
        stress_flag = 0
        stress = []

        for line in lines:
            if stress_flag == 3:
                if line == "\n":
                    stress_flag = 0
                else:
                    stress.append(line.split()[1:4])
            if stress_flag == 2:
                stress_flag = 3
            if stress_flag == 1:
                stress_flag = 2
            if stress_sign in line:
                stress_flag = 1
            if force_start_pattern.match(line):
                force_flag = True
            if force_end_pattern.match(line):
                assert force_flag is True, (
                    force_flag,
                    "there may be errors in this file ",
                )
                force_flag = False
            if force_flag is True:
                force_lines.append(line)
            if energy_pattern_1.match(line):
                energy = (
                    float(energy_pattern_1.match(line).groupdict()["number"]) * AU_TO_EV
                )
                # print('1to', energy)
            if energy_pattern_2.match(line):
                energy = (
                    float(energy_pattern_2.match(line).groupdict()["number"]) * AU_TO_EV
                )
            if cell_length_pattern.match(line):
                cell_A = (
                    float(cell_length_pattern.match(line).groupdict()["A"]) * AU_TO_ANG
                )
                cell_B = (
                    float(cell_length_pattern.match(line).groupdict()["B"]) * AU_TO_ANG
                )
                cell_C = (
                    float(cell_length_pattern.match(line).groupdict()["C"]) * AU_TO_ANG
                )
                cell_flag += 1
            if cell_angle_pattern.match(line):
                cell_alpha = np.deg2rad(
                    float(cell_angle_pattern.match(line).groupdict()["alpha"])
                )
                cell_beta = np.deg2rad(
                    float(cell_angle_pattern.match(line).groupdict()["beta"])
                )
                cell_gamma = np.deg2rad(
                    float(cell_angle_pattern.match(line).groupdict()["gamma"])
                )
                cell_flag += 1
            if print_level_pattern.match(line):
                print_level = print_level_pattern.match(line).groupdict()["print_level"]
                print_level_flag += 1
            if cell_a_pattern.match(line):
                cell_ax = float(cell_a_pattern.match(line).groupdict()["ax"])
                cell_ay = float(cell_a_pattern.match(line).groupdict()["ay"])
                cell_az = float(cell_a_pattern.match(line).groupdict()["az"])
                cell_flag += 1
            if cell_b_pattern.match(line):
                cell_bx = float(cell_b_pattern.match(line).groupdict()["bx"])
                cell_by = float(cell_b_pattern.match(line).groupdict()["by"])
                cell_bz = float(cell_b_pattern.match(line).groupdict()["bz"])
                cell_flag += 1
            if cell_c_pattern.match(line):
                cell_cx = float(cell_c_pattern.match(line).groupdict()["cx"])
                cell_cy = float(cell_c_pattern.match(line).groupdict()["cy"])
                cell_cz = float(cell_c_pattern.match(line).groupdict()["cz"])
                cell_flag += 1

            if atomic_kinds_pattern.match(line):
                akind = atomic_kinds_pattern.match(line).groupdict()["akind"]
                atomic_kinds.append(akind)
        if print_level_flag == 1:
            self.print_level = print_level
            if print_level == "LOW":
                raise RuntimeError(
                    "please provide cp2k output with higher print level(at least MEDIUM)"
                )

        if cell_flag == 2:
            self.cell = cell_to_low_triangle(
                cell_A, cell_B, cell_C, cell_alpha, cell_beta, cell_gamma
            )
        elif cell_flag == 5:
            self.cell = np.asarray(
                [
                    [cell_ax, cell_ay, cell_az],
                    [cell_bx, cell_by, cell_bz],
                    [cell_cx, cell_cy, cell_cz],
                ]
            ).astype("float64")
        if atomic_kinds:
            self.atomic_kinds = atomic_kinds
        # print(self.atomic_kinds)
        # lx = cell_A
        # xy = cell_B * np.cos(cell_gamma)
        # xz = cell_C * np.cos(cell_beta)
        # ly = cell_B* np.sin(cell_gamma)
        # yz = (cell_B*cell_C*np.cos(cell_alpha)-xy*xz)/ly
        # lz = np.sqrt(cell_C**2-xz**2-yz**2)
        # self.cell = [[lx, 0 , 0],
        #         [xy, ly, 0 ],
        #         [xz, yz, lz]]

        element_index = -1
        element_dict = OrderedDict()
        atom_types_idx_list = []
        forces_list = []
        for line in force_lines[3:]:
            line_list = line.split()
            # print(line_list)
            if element_dict.get(line_list[1]):
                element_dict[line_list[1]][1] += 1
            else:
                element_index += 1
                element_dict[line_list[1]] = [element_index, 1]
            atom_types_idx_list.append(element_dict[line_list[1]][0])
            forces_list.append(
                [
                    float(line_list[3]) * AU_TO_EV_EVERY_ANG,
                    float(line_list[4]) * AU_TO_EV_EVERY_ANG,
                    float(line_list[5]) * AU_TO_EV_EVERY_ANG,
                ]
            )
        # print(atom_types_idx_list)
        # atom_names=list(element_dict.keys())
        atom_names = self.atomic_kinds
        atom_numbs = []

        GPa = PressureConversion("eV/angstrom^3", "GPa").value()
        if stress:
            stress = np.array(stress)
            stress = stress.astype("float64")
            stress = stress[np.newaxis, :, :]
            # stress to virial conversion, default unit in cp2k is GPa
            # note the stress is virial = stress * volume
            virial = stress * np.linalg.det(self.cell) / GPa
            virial = virial.squeeze()
        else:
            virial = None
        for ii in element_dict.keys():
            atom_numbs.append(element_dict[ii][1])
        # print(atom_numbs)
        info_dict["atom_names"] = atom_names
        info_dict["atom_numbs"] = atom_numbs
        info_dict["atom_types"] = np.asarray(atom_types_idx_list)
        info_dict["print_level"] = self.print_level
        info_dict["cells"] = np.asarray([self.cell]).astype("float64")
        info_dict["energies"] = np.asarray([energy]).astype("float64")
        info_dict["forces"] = np.asarray([forces_list]).astype("float64")
        if virial is not None:
            info_dict["virials"] = np.asarray([virial]).astype("float64")
        return info_dict

    def handle_single_xyz_frame(self, lines):
        info_dict = {}
        atom_num = int(lines[0].strip("\n").strip())
        if len(lines) != atom_num + 2:
            raise RuntimeError(
                f"format error, atom_num=={atom_num}, {len(lines)}!=atom_num+2"
            )
        data_format_line = lines[1].strip("\n").strip() + " "
        prop_pattern = re.compile(r"(?P<prop>\w+)\s*=\s*(?P<number>.*?)[, ]")
        prop_dict = dict(prop_pattern.findall(data_format_line))

        energy = 0
        if prop_dict.get("E"):
            energy = float(prop_dict.get("E")) * AU_TO_EV
            # info_dict['energies'] = np.array([prop_dict['E']]).astype('float64')

        element_index = -1
        element_dict = OrderedDict()
        atom_types_list = []
        coords_list = []
        for line in lines[2:]:
            line_list = line.split()
            if element_dict.get(line_list[0]):
                element_dict[line_list[0]][1] += 1
            else:
                element_index += 1
                element_dict[line_list[0]] = [element_index, 1]
            atom_types_list.append(element_dict[line_list[0]][0])
            # coords_list.append([float(line_list[1])*AU_TO_ANG,
            #     float(line_list[2])*AU_TO_ANG,
            #     float(line_list[3])*AU_TO_ANG])
            coords_list.append(
                [float(line_list[1]), float(line_list[2]), float(line_list[3])]
            )
        atom_names = list(element_dict.keys())
        atom_numbs = []
        for ii in atom_names:
            atom_numbs.append(element_dict[ii][1])
        # info_dict['atom_names'] = atom_names
        # info_dict['atom_numbs'] = atom_numbs
        # info_dict['atom_types'] = np.asarray(atom_types_list)
        info_dict["coords"] = np.asarray([coords_list]).astype("float64")
        info_dict["energies"] = np.array([energy]).astype("float64")
        info_dict["orig"] = np.zeros(3)
        return info_dict


# %%


def get_frames(fname):
    coord_flag = False
    force_flag = False
    stress_flag = False
    eV = EnergyConversion("hartree", "eV").value()
    angstrom = LengthConversion("bohr", "angstrom").value()
    GPa = PressureConversion("eV/angstrom^3", "GPa").value()
    atom_symbol_idx_list = []
    atom_symbol_list = []
    cell = []
    coord = []
    force = []
    stress = []

    fp = open(fname)
    # check if output is converged, if not, return sys = 0
    content = fp.read()
    count = content.count("SCF run converged")
    if count == 0:
        fp.close()
        return [], [], [], [], [], [], [], None

    # search duplicated header
    fp.seek(0)
    header_idx = []
    for idx, ii in enumerate(fp):
        if "Multiplication driver" in ii:
            header_idx.append(idx)

    # parse from last header
    fp.seek(0)
    for idx, ii in enumerate(fp):
        if idx > header_idx[-1]:
            if "CELL| Vector" in ii:
                cell.append(ii.split()[4:7])
            if "Atomic kind:" in ii:
                atom_symbol_list.append(ii.split()[3])

            # beginning of coords block
            if "Atom  Kind  Element" in ii or "Atom Kind Element" in ii:
                coord_flag = True
            # parse coords lines
            elif coord_flag:
                if ii == "\n":
                    coord_flag = len(coord) == 0  # skip empty line at the beginning
                else:
                    coord.append(ii.split()[4:7])
                    atom_symbol_idx_list.append(ii.split()[1])

            if "ENERGY|" in ii:
                energy = ii.split()[8]
            if " Atom   Kind " in ii:
                force_flag = True
                force_idx = idx
            if force_flag:
                if idx > force_idx:
                    if "SUM OF ATOMIC FORCES" in ii:
                        force_flag = False
                    else:
                        force.append(ii.split()[3:6])
            # add reading stress tensor
            if "STRESS TENSOR [GPa" in ii:
                stress_flag = True
                stress_idx = idx
            if stress_flag:
                if idx > stress_idx + 2:
                    if ii == "\n":
                        stress_flag = False
                    else:
                        stress.append(ii.split()[1:4])

    fp.close()
    assert coord, "cannot find coords"
    assert energy, "cannot find energies"
    assert force, "cannot find forces"

    # conver to float array and add extra dimension for nframes
    cell = np.array(cell)
    cell = cell.astype("float64")
    cell = cell[np.newaxis, :, :]
    coord = np.array(coord)
    coord = coord.astype("float64")
    coord = coord[np.newaxis, :, :]
    atom_symbol_idx_list = np.array(atom_symbol_idx_list)
    atom_symbol_idx_list = atom_symbol_idx_list.astype(int)
    atom_symbol_idx_list = atom_symbol_idx_list - 1
    atom_symbol_list = np.array(atom_symbol_list)
    atom_symbol_list = atom_symbol_list[atom_symbol_idx_list]
    force = np.array(force)
    force = force.astype("float64")
    force = force[np.newaxis, :, :]

    # virial is not necessary
    if stress:
        stress = np.array(stress)
        stress = stress.astype("float64")
        stress = stress[np.newaxis, :, :]
        # stress to virial conversion, default unit in cp2k is GPa
        # note the stress is virial = stress * volume
        virial = stress * np.linalg.det(cell[0]) / GPa
    else:
        virial = None

    # force unit conversion, default unit in cp2k is hartree/bohr
    force = force * eV / angstrom
    # energy unit conversion, default unit in cp2k is hartree
    energy = float(energy) * eV
    energy = np.array(energy).astype("float64")
    energy = energy[np.newaxis]

    tmp_names, symbol_idx = np.unique(atom_symbol_list, return_index=True)
    atom_types = []
    atom_numbs = []
    # preserve the atom_name order
    atom_names = atom_symbol_list[np.sort(symbol_idx, kind="stable")]
    for jj in atom_symbol_list:
        for idx, ii in enumerate(atom_names):
            if jj == ii:
                atom_types.append(idx)
    for idx in range(len(atom_names)):
        atom_numbs.append(atom_types.count(idx))

    atom_types = np.array(atom_types)

    return list(atom_names), atom_numbs, atom_types, cell, coord, energy, force, virial


# %%
