import xml.etree.ElementTree as ET
import os
import argparse, re
import numpy as np
import copy
import Molcas_Methods
from sys import argv

script, first, second, third = argv
cml_file = first
scratchdir = third
global name
name = second

def Get_Order(cml):
    cml = cml.replace(".cml","")
    cml = cml.split("_")
    return int(cml[2])


def Get_Geom_String(cml):
    tree = ET.parse(cml)
    root = tree.getroot()

    geom_string = ''
    for atom in root[0]:
        if 'X' not in str(atom.attrib['elementType']):
            geom_string = geom_string+" "+str(atom.attrib['elementType'])+" "+str(atom.attrib['x3'])+" "+str(atom.attrib['y3'])+" "+str(atom.attrib['z3'])+"\n"
    geom_string+='units angstrom'
    return geom_string

def G_Thread(cml_file):
    np.set_printoptions(threshold = 10000000000000)
    tree = ET.parse(cml_file)
    root = tree.getroot()
    method = root.find('theory').attrib['name']
    basis = root.find('basis').attrib['name']
    sign = int(root.find('sign').attrib['name'])
    chargemult = str(root.find('chargemult').attrib['cm'])
    charge = int(chargemult.split(' ')[1])
    import psi4
    order = Get_Order(cml_file)
    geom = Get_Geom_String(cml_file)
    if method=='mcscf' or method=='MCSCF' or method=='casscf' or method=='CASSCF':
        grad, energy = Molcas_Methods.Molcas_Eval(cml_file, scratchdir)
    else:
        psi4.geometry(geom)
        psi4.core.be_quiet()
        psi4.set_memory('1 GB')
        psi4.core.be_quiet()
        grad, wfn = psi4.gradient(method+'/'+basis, return_wfn=True, CFOUR_CHARGE=charge)
        energy = wfn.energy()*order
    grad = np.asarray(grad)
    true_grad = Interpret(grad, cml_file)*order
    return tuple([sign*true_grad, sign*energy])


def Interpret(raw_grad, cml_file):
    tree = ET.parse(name)
    root = tree.getroot()
    true_grad = []
    for atom in root[0]:
        true_grad.append([0, 0, 0])
    tree2 = ET.parse(cml_file)
    root2 = tree2.getroot()
    index = 0
    for atom in range(0, len(root2[0])):
        if 'X' not in str(root2[0][atom].attrib['elementType']) and atom<len(root[0]):
            true_grad[atom][0] = (raw_grad[index][0])
            true_grad[atom][1] = (raw_grad[index][1])
            true_grad[atom][2] = (raw_grad[index][2])
            index = index+1
        elif 'X' not in str(root2[0][atom].attrib['elementType']):
            info_atom = atom+1
            line = str(root2[0][info_atom].attrib['elementType']).split('_')
            stay = int(line[1])
            leave = int(line[3])
            factor = float(line[4])
            #print("Factor = "+str(factor))
            true_grad[stay][0] += (raw_grad[index][0])*(1-factor)
            true_grad[stay][1] += (raw_grad[index][1])*(1-factor)
            true_grad[stay][2] += (raw_grad[index][2])*(1-factor)
            true_grad[leave][0] += (raw_grad[index][0])*factor
            true_grad[leave][1] += (raw_grad[index][1])*factor
            true_grad[leave][2] += (raw_grad[index][2])*factor
            index = index+1

    return (np.asarray(true_grad))

calc = G_Thread(cml_file)
numero = cml_file.split('_')[1]
results = open(scratchdir+'/'+'res/'+str(numero), 'w')
results.write(str(calc[1])+'\n')
results.write(str(calc[0]))
