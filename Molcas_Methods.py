import os
import xml.etree.ElementTree as ET
import h5py

#Take a cml file and get a CASSCF energy/gradient
def Molcas_Eval(cml, scratchdir):
    tree = ET.parse(cml)
    root = tree.getroot()
    cml = cml.split('/')[-1]
    Make_Subdir(cml, scratchdir)
    Write_XYZ(cml, root)
    Write_Input(cml, root)
    cml_short = cml.replace('.cml', '')
    os.system('pymolcas '+cml_short+'.input')
    gradh5 = h5py.File(cml_short+'.slapaf.h5', 'r')
    gradient = -gradh5['FORCES'][0]
    energy = gradh5['ENERGIES'][0]
    os.chdir('..')
    os.system('rm -r -f '+str(cml_short))
    return ([gradient, energy])
    
#Construct a sub-directory and move into it.
def Make_Subdir(cml, scratchdir):
    cml_short = cml.replace('.cml', '')
    if os.path.exists(cml_short):
        os.system('rm -r -f '+str(cml_short))
    os.system('mkdir '+str(cml_short))
    os.chdir(cml_short)

#Write an XYZ file to hit with MOLCAS
def Write_XYZ(cml, root):
    cml_short = cml.replace('.cml', '')
    os.system('pwd')
    xyz = open(cml_short+'.xyz', 'w')
    n = 0
    for atom in root[0]:
        if 'X' not in str(atom.attrib['elementType']):
            n+=1
    xyz.write(str(n)+'\n')
    xyz.write('Angstroms\n')
    for atom in root[0]:
        if 'X' not in str(atom.attrib['elementType']):
            xyz.write(str(atom.attrib['elementType'])+' '+str(atom.attrib['x3'])+' '+str(atom.attrib['y3'])+' '+str(atom.attrib['z3'])+'\n')
    xyz.close()
            
def Write_Input(cml, root):
    cml_short = cml.replace('.cml', '')
    basis = root.find('basis').attrib['name']
    chargemult = str(root.find('chargemult').attrib['cm'])
    charge = chargemult.split(' ')[0]
    mult = chargemult.split(' ')[1]
    inp = open(cml_short+'.input', 'w')
    inp.write('&GATEWAY\n')
    inp.write(' coord='+cml_short+'.xyz\n')
    inp.write(' basis='+str(basis)+'\n')
    inp.write('&SEWARD\n')
    inp.write('&CASSCF\n')
    inp.write(' CHARGE = '+str(charge)+'\n')
    inp.write(' SPIN = '+str(mult)+'\n')
    inp.write('&ALASKA\n')
    inp.write('&SLAPAF\n') 
    
