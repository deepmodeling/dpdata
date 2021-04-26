def to_gjf_string(system, frame_idx=0, header="", title="", foot="", charge=0, mult=1):
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