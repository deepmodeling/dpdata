import numpy as np
import re

latt_patt="\|\s+([0-9]{1,}[.][0-9]*)\s+([0-9]{1,}[.][0-9]*)\s+([0-9]{1,}[.][0-9]*)"
pos_patt_first="\|\s+[0-9]{1,}[:]\s\w+\s(\w+)(\s.*[-]?[0-9]{1,}[.][0-9]*)(\s+[-]?[0-9]{1,}[.][0-9]*)(\s+[-]?[0-9]{1,}[.][0-9]*)"
pos_patt_other="\s+[a][t][o][m]\s+([-]?[0-9]{1,}[.][0-9]*)\s+([-]?[0-9]{1,}[.][0-9]*)\s+([-]?[0-9]{1,}[.][0-9]*)\s+(\w{1,2})"
force_patt="\|\s+[0-9]{1,}\s+([-]?[0-9]{1,}[.][0-9]*[E][+-][0-9]{1,})\s+([-]?[0-9]{1,}[.][0-9]*[E][+-][0-9]{1,})\s+([-]?[0-9]{1,}[.][0-9]*[E][+-][0-9]{1,})"
eng_patt="Total energy uncorrected.*([-]?[0-9]{1,}[.][0-9]*[E][+-][0-9]{1,})\s+eV"
#atom_numb_patt="Number of atoms.*([0-9]{1,})"

debug = False
def get_info (lines, type_idx_zero = False) :

    atom_types = []
    atom_names = []
    cell = []
    atom_numbs = None
    _atom_names = []

    contents="\n".join(lines)
    #cell
    #_tmp=re.findall(latt_patt,contents)
    #for ii in _tmp:
    #    vect=[float(kk) for kk in ii]
    #    cell.append(vect)
    #------------------
    for ln,l in enumerate(lines):
        if l.startswith('  | Unit cell'):
            break
    _tmp=lines[ln+1:ln+4]
    for ii in _tmp:
        v_str=ii.split('|')[1].split()
        vect=[float(kk) for kk in v_str]
        cell.append(vect)

    _tmp=re.findall(pos_patt_first,contents)
    for ii in _tmp:
        _atom_names.append(ii[0])
    atom_names=[]
    for ii in _atom_names:
        if not ii in atom_names:
           atom_names.append(ii)
    
    atom_numbs =[_atom_names.count(ii) for ii in atom_names] 
    if type_idx_zero :
       type_map=dict(zip(atom_names,range(len(atom_names)))) 
    else:
       type_map=dict(zip(atom_names,range(1,len(atom_names)+1))) 
    atom_types=list(map(lambda  k: type_map[k], _atom_names)) 
    assert(atom_numbs is not None), "cannot find ion type info in aims output"
     

    return [cell, atom_numbs, atom_names, atom_types ]


def get_fhi_aims_block(fp) :
    blk = []
    for ii in fp :
        if not ii :
            return blk
        blk.append(ii.rstrip('\n'))
        if 'Begin self-consistency loop: Re-initialization' in ii:
            return blk
    return blk

def get_frames (fname, md=True, begin = 0, step = 1) :
    fp = open(fname)
    blk = get_fhi_aims_block(fp)
    ret = get_info(blk, type_idx_zero = True)

    cell, atom_numbs, atom_names, atom_types =ret[0],ret[1],ret[2],ret[3]
    ntot = sum(atom_numbs)

    all_coords = []
    all_cells = []
    all_energies = []
    all_forces = []
    all_virials = []    

    cc = 0
    while len(blk) > 0 :
        if debug:
           with open(str(cc),'w') as f:
                f.write('\n'.join(blk))
        if cc >= begin and (cc - begin) % step == 0 :
            if cc==0:
                coord, _cell, energy, force, virial, is_converge = analyze_block(blk, first_blk=True, md=md)
            else:
                coord, _cell, energy, force, virial, is_converge = analyze_block(blk, first_blk=False)
            if is_converge : 
                if len(coord) == 0:
                    break
                all_coords.append(coord)

                if _cell:
                   all_cells.append(_cell)
                else:
                   all_cells.append(cell)

                all_energies.append(energy)
                all_forces.append(force)
                if virial is not None :
                    all_virials.append(virial)
        blk = get_fhi_aims_block(fp)
        cc += 1
        
    if len(all_virials) == 0 :
        all_virials = None
    else :
        all_virials = np.array(all_virials)
    fp.close()
    return atom_names, atom_numbs, np.array(atom_types), np.array(all_cells), np.array(all_coords), np.array(all_energies), np.array(all_forces), all_virials


def analyze_block(lines, first_blk=False, md=True) :
    coord = []
    cell = []
    energy = None
    force = []
    virial = None
    atom_names=[]
    _atom_names=[]

    contents="\n".join(lines)
    try:
       natom=int(re.findall("Number of atoms.*([0-9]{1,})",lines)[0])
    except:
       natom=0

    if first_blk:

       if md:
          _tmp=re.findall(pos_patt_other,contents)[:]
          for ii in _tmp[slice(int(len(_tmp)/2),len(_tmp))]:
              coord.append([float(kk) for kk in ii[:-1]])
       else:
          _tmp=re.findall(pos_patt_first,contents)
          for ii in _tmp:
              coord.append([float(kk) for kk in ii[1:]])
    else:
       _tmp=re.findall(pos_patt_other,contents)
       for ii in _tmp:
           coord.append([float(kk) for kk in ii[:-1]])

    _tmp=re.findall(force_patt,contents)
    for ii in _tmp:
        force.append([float(kk) for kk in ii])

    if "Self-consistency cycle converged" in contents:
       is_converge=True
    else:
       is_converge=False

    try:
      _eng_patt=re.compile(eng_patt)
      energy=float(_eng_patt.search(contents).group().split()[-2])
    except:
     energy=None
    
    if not energy:
       is_converge = False

    if energy:
       assert((force is not None) and len(coord) > 0 )

    return coord, cell, energy, force, virial, is_converge

if __name__=='__main__':
  import sys
  ret=get_frames (sys.argv[1], begin = 0, step = 1)
  print(ret)
