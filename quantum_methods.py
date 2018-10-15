import xml.etree.ElementTree as ET
import os
import argparse, re
import numpy as np
import copy
import multiprocessing
from multiprocessing import Pool


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
    return geom_string

def Compute_Gradient(args, root):
    gradient = []
    energy = 0
    name = str(args['og_name']).replace('[','')
    name = name.replace(']','')
    name = str(args['og_name']).replace('[','')
    name = name.replace(']','')
    name = name.replace('\'','')
    tree = ET.parse(name)
    root = tree.getroot()
    for atom in root[0]:
        gradient.append([0.0,0.0,0.0])
    gradient = np.asarray(gradient)
    if os.path.exists(args['scratch']+'/res'):
        os.system('rm -r '+str(args['scratch']+'/res'))
    os.system('mkdir '+str(args['scratch']+'/res'))
    for cml_file in os.listdir(args['scratch']+"/cmls"):
        if ".cml" in cml_file:
            cml_file_name = args['scratch']+"/cmls/"+cml_file
            os.system('python Grad_Standalone.py '+str(cml_file_name)+' '+str(name)+' '+args['scratch'])
    for file in os.listdir(args['scratch']+"/res"):
        frag = open(args['scratch']+"/res/"+file, 'r')
        lines = frag.readlines()
        for line in lines:
            line = line.replace(']','')
            line = line.replace('[','')
            line = line.replace('\n','')
        energy += float(lines[0])
        for line in range(1, gradient.shape[0]+1):
            string = lines[line]
            string = string.replace(']','')
            string = string.replace('[','')
            string = string.replace('\n','')
            string=string.split()
            gradient[line-1][0]+=float(string[0])
            gradient[line-1][1]+=float(string[1])
            gradient[line-1][2]+=float(string[2])
    return tuple([gradient, energy])
