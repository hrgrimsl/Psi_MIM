import argparse,re
import xml.etree.ElementTree as ET
import numpy as np
import cProfile
import time
import os

from frag_classes import *
from frag_methods import *
import cython_pie

def add_up_atoms(frag_list, n_atoms):
    mult = [0]*n_atoms
    for frag in frag_list:
        for atom in frag.atoms:
            mult[atom] += frag.order
    print("Number of times each atom is counted:")
    print(mult)

def Fragment_Verb(args, index):
    file_name = args['name']
    file_name = file_name[0]
    assert(file_name[-3:] == "cml")
    scratchdir = args['scratch']
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
    final_frag_list = cython_pie.cython_pie(primary_frags)
    full_system.make_frag_objects(final_frag_list)
    index += 1
    for frag in full_system.frags:
        if frag.order!=0:
            full_system.append_meta_list(frag, scratchdir, index)
            dst_file = str(scratchdir)+"/"+str(file_name).replace(".cml","")+"_"+str(index)+"_"+str(frag.order)+".cml"
            full_system.write_cml(frag, dst_file, args)
            index += 1
    return index

def Super_Fragment(core_args, dict_list):
    if os.path.exists(core_args['scratch']+'/cmls'):
        os.system(('rm -r '+core_args['scratch']+'/cmls'))
    os.system(('mkdir '+core_args['scratch']+'/cmls'))
    metaname = str(core_args['scratch'])+"/cmls/"+str(core_args['name'][0]).replace('.cml',"_CML_list")
    meta = open(metaname, "w")
    scratchdir = str(core_args['scratch'])
    index = 0
    for dict in dict_list:
        index = Fragment_Verb(dict, index)
