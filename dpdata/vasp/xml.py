#!/usr/bin/env python3

import xml.etree.ElementTree as ET
import numpy as np

def check_name(item, name) :
    assert (item.attrib['name'] == name), "item attrib '%s' dose not math required '%s'" % (item.attrib['name'], name)

def get_varray(varray) :
    array = []
    for vv in varray.findall('v') :
        array.append([ float(ii) for ii in vv.text.split()])
    return np.array(array)

def analyze_atominfo(atominfo_xml) :
    check_name(atominfo_xml.find('array'), 'atoms')
    eles = []
    types = []
    for ii in atominfo_xml.find('array').find('set') :
        eles .append((ii.findall('c')[0].text.strip()))
        types.append(int(ii.findall('c')[1].text))
    uniq_ele = []
    for ii in eles :
        if not(ii  in uniq_ele):
            uniq_ele.append(ii)
    return uniq_ele, types

def analyze_calculation(cc) :
    structure_xml = cc.find('structure')
    check_name(structure_xml.find('crystal').find('varray'), 'basis')
    check_name(structure_xml.find('varray'), 'positions')
    cell = get_varray(structure_xml.find('crystal').find('varray'))
    posi = get_varray(structure_xml.find('varray'))    
    strs = None
    for vv in cc.findall('varray') :
        if vv.attrib['name'] == 'forces' :
            forc = get_varray(vv) 
        elif vv.attrib['name'] == 'stress' :
            strs = get_varray(vv)
    for ii in cc.find('energy').findall('i') :
        if ii.attrib['name'] == 'e_fr_energy' :
            ener = float(ii.text)
    # print(ener)
    # return 'a'
    return posi, cell, ener, forc, strs

def formulate_config(eles, types, posi, cell, ener, forc, strs_) :
    strs = strs_ / 1602
    natoms = len(types)
    ntypes = len(eles)    
    ret = ""
    ret += "#N %d %d\n" % (natoms, ntypes-1)
    ret += "#C "
    for ii in eles :
        ret += ' ' + ii
    ret += '\n'
    ret += "##\n"
    ret += '#X %13.8f %13.8f %13.8f\n' % (cell[0][0], cell[0][1], cell[0][2])
    ret += '#Y %13.8f %13.8f %13.8f\n' % (cell[1][0], cell[1][1], cell[1][2])
    ret += '#Z %13.8f %13.8f %13.8f\n' % (cell[2][0], cell[2][1], cell[2][2])
    ret += "#W 1.0\n"
    ret += "#E %.10f\n" % (ener / natoms)
    ret += '#S %.9e %.9e %.9e %.9e %.9e %.9e\n' % \
           (strs[0][0], strs[1][1], strs[2][2], strs[0][1], strs[1][2], strs[0][2])
    ret += '#F\n'
    for ii in range(natoms) :
        sp = np.matmul(cell.T, posi[ii])
        ret += '%d' % (types[ii]-1)
        ret += ' %12.6f %12.6f %12.6f' % (sp[0], sp[1], sp[2])
        ret += ' %12.6f %12.6f %12.6f' % (forc[ii][0], forc[ii][1], forc[ii][2])
        ret += '\n'            
    return ret

def analyze (fname, type_idx_zero = False, begin = 0, step = 1) :
    """
    can deal with broken xml file
    """
    all_posi = []
    all_cell = []
    all_ener = []
    all_forc = []
    all_strs = []
    cc = 0
    try:
        for event, elem in ET.iterparse(fname):
            if elem.tag == 'atominfo' :
                eles, types = analyze_atominfo(elem)
                types = np.array(types, dtype = int)
                if type_idx_zero :
                    types = types - 1
            if elem.tag == 'calculation' :
                posi, cell, ener, forc, strs = analyze_calculation(elem)
                if cc >= begin and (cc - begin) % step == 0 :
                    all_posi.append(posi)
                    all_cell.append(cell)
                    all_ener.append(ener)
                    all_forc.append(forc)
                    if strs is not None :
                        all_strs.append(strs)                
                cc += 1
    except ET.ParseError:
        return eles, types, np.array(all_cell), np.array(all_posi), np.array(all_ener), np.array(all_forc), np.array(all_strs)
    return eles, types, np.array(all_cell), np.array(all_posi), np.array(all_ener), np.array(all_forc), np.array(all_strs)

