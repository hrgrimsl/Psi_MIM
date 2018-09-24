import argparse,re
import xml.etree.ElementTree as ET
import numpy as np
import cProfile
import time
import os

from frag_classes import *
from frag_methods import *
from cython_pie import *
from cython_pie import cython_pie as cpie


def add_up_atoms(frag_list, n_atoms):
    mult = [0]*n_atoms
    for frag in frag_list:
        for atom in frag.atoms:
            mult[atom] += frag.order
    print(" Number of times each atom is counted:")
    print(mult)

'''
#   Setup input arguments
parser = argparse.ArgumentParser(description='Creates a set of fragments from a cml file',
formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('name', nargs=1, help='CML file to process')
parser.add_argument('-n','--eta', type=int, default=1, help='Which eta to use?', required=False)
parser.add_argument('-r','--radius', type=float, help='Specifies radius for fragmentation (NYI)', required=False)
parser.add_argument('-s','--scratch', type=str, help='Name of scratch directory', default="./scr", required=False)
args = vars(parser.parse_args())
'''
def Fragment(args):
    file_name = args['name']
    file_name = file_name[0]
    assert(file_name[-3:] == "cml")
    print(" Process file:", file_name)

    covrad = form_covalent_radii()
    eta = int(args['eta'])
    full_system = Molecule()
    full_system.parse_cml(file_name)
    full_system.get_prims()
    full_system.get_prim_conns()
    full_system.get_frags(eta)
    full_system.cull_frags()
    full_system.close_frags()
    full_system.cull_frags()
    full_system.finalize_first_frags()
    primary_frags = full_system.construct_frag_dict()

    #Line in question
    final_frag_list = cpie(primary_frags)
    full_system.make_frag_objects(final_frag_list)
    scratchdir = str(args['scratch'])
    if os.path.exists(str(scratchdir)):
        os.system('rm -r '+str(scratchdir))
    os.system('mkdir '+str(scratchdir))
    index = 0
    for frag in full_system.frags:
        if frag.order!=0:
            full_system.append_meta_list(frag, scratchdir, index)
            dst_file = str(scratchdir)+"/"+str(file_name).replace(".cml","")+"_"+str(index)+"_"+str(frag.order)+".cml"
            full_system.write_cml(frag, dst_file)
            index += 1
    print(len(full_system.frags))
