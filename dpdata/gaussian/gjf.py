def to_gjf_string(system, frame_idx=0, header="", title="", foot="", charge=0, mult=1):
    '''
    Generate a string that can be written to a Gaussian input file, please refer to https://gaussian.com/input/?tabid=2 for more information on syntax of Gaussian input files.

    Parameters
    --------
    system : dpdata.System
        The system to write a gaussian input file
    frame_idx : int
        The index of frame to specify molecule geometry
    header : str
        The route section (# lines) and link-0 commands (% commands)
    title : str
        The title section of .gjf files
    foot : str
        The modifications to coordinates, used when Opt=ModRedundant
    charge : int
        The charge of the system
    mult : int
        The multiplicity of the system
    '''
    if not title:
        title = system.formula
    if not header:
        header = "#force B3LYP/6-31G(d)"
    ret = header + "\n\n" + title + "\n\n" + str(charge) + " " + str(mult) + "\n"
    coord = system.data["coords"][frame_idx].reshape(-1, 3)
    for atype, c in zip(system.data["atom_types"], coord):
        ret += system.data["atom_names"][atype] + " "
        ret += " ".join([str(x) for x in c])
        ret += "\n"
    if foot:
        ret += "\n" + foot + "\n"
    else:
        ret += "\n"
    return ret