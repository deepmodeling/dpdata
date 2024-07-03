#!/usr/bin/env python3
# %%
from __future__ import annotations

import re
from collections import OrderedDict

import numpy as np
from ..unit import EnergyConversion, ForceConversion, LengthConversion

e_conv_kcalpermol2eV = EnergyConversion("kcal_mol", "eV").value()
e_conv_au2eV = EnergyConversion("hartree", "eV").value()
f_conv_kcalpermolperang2eVperang = ForceConversion("kcal_mol/angstrom", "eV/angstrom").value()
f_conv_auperang2eVperang = ForceConversion("hartree/angstrom", "eV/angstrom").value()
f_conv_kcalpermolperbohr2eVperang = ForceConversion("kcal_mol/bohr", "eV/angstrom").value()
f_conv_au2eVperang = ForceConversion("hartree/bohr", "eV/angstrom").value()


class QuipGapxyzSystems:
    """deal with QuipGapxyzFile."""

    def __init__(self, file_name):
        self.file_object = open(file_name)
        self.block_generator = self.get_block_generator()

    def __iter__(self):
        return self

    def __next__(self):
        return self.handle_single_xyz_frame(next(self.block_generator))

    def __del__(self):
        self.file_object.close()

    def get_block_generator(self):
        p3 = re.compile(r"^\s*(\d+)\s*")
        while True:
            line = self.file_object.readline()
            if not line:
                break
            if p3.match(line):
                atom_num = int(p3.match(line).group(1))
                lines = []
                lines.append(line)
                for ii in range(atom_num + 1):
                    lines.append(self.file_object.readline())
                if not lines[-1]:
                    raise RuntimeError(
                        "this xyz file may lack of lines, should be {};lines:{}".format(
                            atom_num + 2, lines
                        )
                    )
                yield lines

        
    @staticmethod
    def handle_single_xyz_frame(lines):
        atom_num = int(lines[0].strip("\n").strip())
        if len(lines) != atom_num + 2:
            raise RuntimeError(
                f"format error, atom_num=={atom_num}, {len(lines)}!=atom_num+2"
            )
        data_format_line = lines[1].strip("\n").strip() + " "
        field_value_pattern = re.compile(
            r"(?P<key>\S+)=(?P<quote>[\'\"]?)(?P<value>.*?)(?P=quote)\s+"
        )
        prop_pattern = re.compile(
            r"(?P<key>\w+?):(?P<datatype>[a-zA-Z]):(?P<value>\d+)"
        )

        data_format_list = [
            kv_dict.groupdict()
            for kv_dict in field_value_pattern.finditer(data_format_line)
        ]
        field_dict = {}
        for item in data_format_list:
            field_dict[item["key"]] = item["value"]

        Properties = field_dict["Properties"]
        prop_list = [
            kv_dict.groupdict() for kv_dict in prop_pattern.finditer(Properties)
        ]

        data_lines = []
        for line in lines[2:]:
            data_lines.append(list(filter(bool, line.strip().split())))
        data_array = np.array(data_lines)
        used_colomn = 0

        type_array = None
        coords_array = None
        Z_array = None
        force_array = None
        virials = None
        for kv_dict in prop_list:
            try:
                if kv_dict["key"] == "species":
                    if kv_dict["datatype"] != "S":
                        raise RuntimeError(
                            "datatype for species must be 'S' instead of {}".format(
                                kv_dict["datatype"]
                            )
                        )
                    field_length = int(kv_dict["value"])
                    type_array = data_array[
                        :, used_colomn : used_colomn + field_length
                    ].flatten()
                    used_colomn += field_length
                    continue
                elif kv_dict["key"] == "pos":
                    if kv_dict["datatype"] != "R":
                        raise RuntimeError(
                            "datatype for pos must be 'R' instead of {}".format(
                                kv_dict["datatype"]
                            )
                        )
                    field_length = int(kv_dict["value"])
                    coords_array = data_array[:, used_colomn : used_colomn + field_length]
                    used_colomn += field_length
                    continue
                elif kv_dict["key"] == "Z":
                    if kv_dict["datatype"] != "I":
                        raise RuntimeError(
                            "datatype for pos must be 'R' instead of {}".format(
                                kv_dict["datatype"]
                            )
                        )
                    field_length = int(kv_dict["value"])
                    Z_array = data_array[
                        :, used_colomn : used_colomn + field_length
                    ].flatten()
                    used_colomn += field_length
                    continue
                elif kv_dict["key"] == "tags":
                    if kv_dict["datatype"] != "I":
                        raise RuntimeError(
                            "datatype for tags must be 'I' instead of {}".format(
                                kv_dict["datatype"]
                            )
                        )
                    field_length = int(kv_dict["value"])
                    used_colomn += field_length
                    continue            
                elif kv_dict["key"] == "force" or kv_dict["key"] == "forces" :
                    if kv_dict["datatype"] != "R":
                        raise RuntimeError(
                            "datatype for pos must be 'R' instead of {}".format(
                                kv_dict["datatype"]
                            )
                        )
                    field_length = int(kv_dict["value"])
                    force_array = data_array[:, used_colomn : used_colomn + field_length]
                    used_colomn += field_length
                    continue
            except Exception as e:
                print("unknown field {}".format(kv_dict["key"]), e)
                continue

        type_num_dict = OrderedDict()
        atom_type_list = []
        type_map = {}
        temp_atom_max_index = 0
        if type_array is None:
            raise RuntimeError("type_array can't be None type, check .xyz file")
        for ii in type_array:
            if ii not in type_map:
                type_map[ii] = temp_atom_max_index
                temp_atom_max_index += 1
                temp_atom_index = type_map[ii]
                atom_type_list.append(temp_atom_index)
                type_num_dict[ii] = 1
            else:
                temp_atom_index = type_map[ii]
                atom_type_list.append(temp_atom_index)
                type_num_dict[ii] += 1
        type_num_list = []
        for atom_type, atom_num in type_num_dict.items():
            type_num_list.append((atom_type, atom_num))
        type_num_array = np.array(type_num_list)
        if field_dict.get("virial", None):
            virials = np.array(
                [
                    np.array(
                        list(filter(bool, field_dict["virial"].split(" ")))
                    ).reshape(3, 3)
                ]
            ).astype("float32")
        else:
            virials = None
            
        try:
            e_units = np.array([field_dict["energy-unit"].lower()])
            f_units = np.array([field_dict["force-unit"].lower()])
        except:
            pass
            #print('No units information contained.')


        info_dict = {}
        info_dict["nopbc"] = False
        info_dict["atom_names"] = list(type_num_array[:, 0])
        info_dict["atom_numbs"] = list(type_num_array[:, 1].astype(int))
        info_dict["atom_types"] = np.array(atom_type_list).astype(int)

        info_dict["coords"] = np.array([coords_array]).astype("float32")
        '''
        info_dict["energies"] = np.array([field_dict["energy"]]).astype("float32")
        info_dict["forces"] = np.array([force_array]).astype("float32")
        '''

        try:
            if e_units == "kcal/mol":
                info_dict["energies"] = np.array([field_dict["energy"]]).astype("float32") * e_conv_kcalpermol2eV
            elif e_units in ["hartree", "au", "a.u."]:
                info_dict["energies"] = np.array([field_dict["energy"]]).astype("float32") * e_conv_au2eV
            elif e_units == "ev":
                info_dict["energies"] = np.array([field_dict["energy"]]).astype("float32")
            else: info_dict["energies"] = np.array([field_dict["energy"]]).astype("float32")
        except Exception:
            try:
                possible_fields = [
                    "energy",
                    "energies",
                    "Energies",
                    "potential-energy.energy",
                    "Energy"
                ]
                for key in possible_fields:
                    if key in field_dict:
                        info_dict["energies"] = np.array([field_dict[key]], dtype="float32")
                        break
                else:
                    raise ValueError("No valid energy field found in field_dict.")
            except KeyError:
                raise ValueError("Error while accessing energy fields in field_dict.")
        
        try:
            if f_units == "kcal/mol/angstrom":
                info_dict["forces"] = np.array([force_array]).astype("float32") * f_conv_kcalpermolperang2eVperang
            elif f_units == "hartree/angstrom" or f_units == "hartree/ang" or f_units == "hartree/ang.":
                info_dict["forces"] = np.array([force_array]).astype("float32") * f_conv_auperang2eVperang
            elif f_units == "kcal/mol/bohr":    
                info_dict["forces"] = np.array([force_array]).astype("float32") * f_conv_kcalpermolperbohr2eVperang
            elif f_units == "kcal/mol/bohr":    
                info_dict["forces"] = np.array([force_array]).astype("float32") * f_conv_au2eVperang
            elif f_units == "ev/angstrom" or f_units == "ev/ang" or f_units == "ev/ang.":
                info_dict["forces"] = np.array([force_array]).astype("float32") 
            else: info_dict["forces"] = np.array([force_array]).astype("float32")
        except:
            info_dict["forces"] = np.array([force_array]).astype("float32") 

        if virials is not None:
            info_dict["virials"] = virials
        info_dict["orig"] = np.zeros(3)
        if "Lattice" in field_dict and field_dict["Lattice"].strip():
            lattice_values = list(filter(bool, field_dict["Lattice"].split(" ")))
            info_dict["cells"] = np.array(
                [np.array(lattice_values).reshape(3, 3)]
            ).astype("float32")
        else:
            lattice_values = np.array([[100.0, 0.0, 0.0], [0.0, 100.0, 0.0], [0.0, 0.0, 100.0]])

            info_dict["nopbc"] = True
            info_dict["cells"] = np.array(
                [np.array(lattice_values).reshape(3, 3)]
            ).astype("float32")

        return info_dict
