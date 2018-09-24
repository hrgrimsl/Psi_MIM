import xml.etree.ElementTree as ET
import os
import argparse, re
import numpy as np
import copy
import multiprocessing
from multiprocessing import Pool

'''
TLDR:

Call

quantum_methods.Compute_Energy(args)

and

quantum_methods.Compute_Gradient(args)

after fragmenting a molecule with frag_script from your wrapper or whatever

args looks like the library below

Gradient will be given as a numpy array that should play nice with psi4

'''
#Argument library for standalone running.
'''
args = {
'basis' : 'STO-3G',
'method' : 'scf',
'name' : 'cluster.cml',
'scratch': './scr'
}
'''


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

def E_Thread(cml_file):
    import psi4
    order = Get_Order(cml_file)
    geom = Get_Geom_String(cml_file)
    psi4.set_options({'FAIL_ON_MAXITER': False})
    psi4.geometry(geom)
    psi4.core.be_quiet()
    if 'bad' in cml_file or 'big' in cml_file:
        method = args1['theoryw']
    if 'good' in cml_file:
        method = args1['theorym']
    basis = args1['basis']
    if 'smallbad' in cml_file:
        return -psi4.energy(method+"/"+basis)*order
    else:
        return psi4.energy(method+"/"+basis)*order
    #energy2 = -psi4.energy(args1['theorym']+"/"+args1['basis'],molecule = dict[cml_file])+psi4.energy(args['theoryw']+"/"+args['basis'],molecule = dict[cml_file])
    #efile = open(args['scratch']+"/efile", 'a')
    #efile.write(str(energy2*order)+"\n")


def Compute_Energy(args):
    global args1
    args1 = args
    procs = []
    dict = {}
    jobs = []
    efile = open(args['scratch']+"/efile",'w')
    pool = Pool(8)
    for cml_file in os.listdir(args['scratch']+"/big"):
        if ".cml" in cml_file:
            cml_file = args['scratch']+"/big/"+cml_file
            jobs.append(cml_file)
            #gradient = gradient + G_Thread(args, cml_file, root)
    for cml_file in os.listdir(args1['scratch']+"/smallbad"):
        if ".cml" in cml_file:
            cml_file = args['scratch']+"/smallbad/"+cml_file
            jobs.append(cml_file)
    for cml_file in os.listdir(args1['scratch']+"/smallgood"):
        if ".cml" in cml_file:
            cml_file = args['scratch']+"/smallgood/"+cml_file
            jobs.append(cml_file)
    results = pool.map(E_Thread, jobs)
    pool.close()
    pool.join()
    energy = 0.0
    for result in results:
        energy = energy+result

    '''
    for cml_file in os.listdir(args['scratch']):
        if ".cml" in cml_file:
            proc = multiprocessing.Process(target=E_Thread, args=(dict, args, cml_file))
            procs.append(proc)
            proc.start()
    for proc in procs:
        proc.join()

    energy = 0.0
    efile = open(args['scratch']+"/efile",'r')
    for e_i in efile.readlines():
        e_i = e_i.replace("\n", "")
        energy = energy+float(e_i)


    for cml_file in os.listdir(args['scratch']):
        if ".cml" in cml_file:
            order = Get_Order(args['scratch']+"/"+cml_file)
            geom = Get_Geom_String(args['scratch']+"/"+cml_file)
            psi4.set_options({'FAIL_ON_MAXITER': False})
            psi4.geometry(geom)
            energy2 = psi4.energy(args['method']+"/"+args['basis'])
            energy = float(energy2)*order+energy
    '''
    return energy

def Interpret(raw_grad, cml_file, args, root):
    true_grad = []
    for atom in root[0]:
        true_grad.append([0, 0, 0])
    tree2 = ET.parse(cml_file)
    root2 = tree2.getroot()
    index = 0
    #print("raw grad")
    #print(raw_grad)
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
            #print("Interpreted = ")
            #print(np.asarray(true_grad))
    return (np.asarray(true_grad))

def G_Thread(cml_file):
    import psi4
    order = Get_Order(cml_file)
    geom = Get_Geom_String(cml_file)
    psi4.set_options({'FAIL_ON_MAXITER': False})
    psi4.geometry(geom)
    psi4.core.be_quiet()
    if 'bad' in str(cml_file) or 'big' in cml_file:
        method = args1['theoryw']
    if 'good' in str(cml_file):
        method = args1['theorym']
    grad, wfn = psi4.gradient(method+'/'+args1['basis'], return_wfn=True)
    energy = wfn.energy()*order
    grad = np.asarray(grad)
    #print(wfn.get_variable('CURRENT ENERGY'))
    #grad = np.asarray(psi4.gradient(method+"/"+args1['basis']))
    #print(energy)
    #gfile = open(args['scratch']+"/gfile", 'a')
    #gfile.write(str(energy2*order)+"\n")
    true_grad = Interpret(grad, cml_file, args1, root1)*order
    if 'smallbad' in cml_file:
        return tuple([-1*true_grad, -1*energy])
    return tuple([true_grad, energy])

def Compute_Gradient(args, root):
    global root1
    root1 = root
    global args1
    args1 = args
    gradient = []
    energy = 0
    name = str(args1['name']).replace('[','')
    name = name.replace(']','')
    name = str(args1['name']).replace('[','')
    name = name.replace(']','')
    name = name.replace('\'','')
    tree = ET.parse(name)
    root = tree.getroot()
    for atom in root[0]:
        gradient.append([0.0,0.0,0.0])
    gradient = np.asarray(gradient)

    pool = Pool(8)
    i = 0
    jobs = []
    for cml_file in os.listdir(args['scratch']+"/big"):
        if ".cml" in cml_file:
            cml_file = args['scratch']+"/big/"+cml_file
            jobs.append(cml_file)
            #gradient = gradient + G_Thread(args, cml_file, root)
    for cml_file in os.listdir(args1['scratch']+"/smallbad"):
        if ".cml" in cml_file:
            cml_file = args['scratch']+"/smallbad/"+cml_file
            jobs.append(cml_file)
    for cml_file in os.listdir(args1['scratch']+"/smallgood"):
        if ".cml" in cml_file:
            cml_file = args['scratch']+"/smallgood/"+cml_file
            jobs.append(cml_file)
    results = pool.map(G_Thread, jobs)
    pool.close()
    pool.join()
    '''
    for cml_file in os.listdir(args['scratch']):
        if ".cml" in cml_file:
            gradient = gradient + G_Thread(cml_file)
    '''
    for result in results:
        gradient += result[0]
        energy += result[1]
    return tuple([gradient, energy])
