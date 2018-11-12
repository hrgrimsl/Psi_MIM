import xml.etree.ElementTree as ET
import psi4

def Get_Geom_String(cml):
    tree = ET.parse(cml)
    root = tree.getroot()
    geom_string = """"""
    cm = root.find('chargemult').attrib['cm'].split()
    for atom in root[0]:
        if 'X' not in str(atom.attrib['elementType']):
            geom_string = geom_string+str(atom.attrib['elementType'])+" "+str(atom.attrib['x3'])+" "+str(atom.attrib['y3'])+" "+str(atom.attrib['z3'])+"\n"
    geom_string+=cm[0]+' '
    geom_string+=cm[1]+'\n'
    geom_string+="symmetry c1"
    return (geom_string, cm[1])

def Parse(cml_file):
    geom, mult = Get_Geom_String(cml_file)
    print(mult)
    print(geom)
    psi4.geometry(geom)
    psi4.set_options({'reference': 'rohf'})
    eng, wfn = psi4.energy('scf/sto-3g', return_wfn=True)
    basis_ = wfn.basisset()
    C_occ = wfn.Ca_subset("AO", "OCC")
    total = len(psi4.core.Matrix.to_array(C_occ)[0])
    active = int(mult)-1
    docc = total-active
    return (docc, active, geom)
