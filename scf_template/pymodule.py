#
# @BEGIN LICENSE
#
# scf_template by Psi4 Developer, a plugin to:
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import psi4
import os
import numpy as np
import time
import psi4.driver.p4util as p4util
from psi4.driver.procrouting import proc_util
import quantum_methods
import frag_script
import argparse, re
import xml.etree.ElementTree as ET





def run_scf_template(name, **kwargs):
    return psi4.core.get_variable('CURRENT ENERGY')

def run_scf_template_grad(name, **kwargs):
    np.set_printoptions(precision=11)
    os.system('cd ..')
    si_file = './SI.txt'
    sifile = open(str(si_file), 'r')
    lines = sifile.readlines()
    for i in range(0, len(lines)):
        lines[i] = (str(lines[i])).strip('\n')
    system_name = str(lines[0])
    file_name = str(system_name) + ".cml"
    etas = lines[1].split("/")
    lil_eta = str(etas[0])
    big_eta = str(etas[1])
    mets = lines[2].split("/")
    lil_met = str(mets[0])
    big_met = str(mets[1])
    bases = lines[3].split("/")
    lil_bas = str(bases[0])
    big_bas = str(bases[1])
    scratchdir = str(lines[4])

    args = {
    'basis': lil_bas,
    'name': [file_name],
    'etam': lil_eta,
    'etaw': big_eta,
    'theorym': lil_met,
    'theoryw': big_met,
    'scratch': scratchdir,
    'og_name': file_name
    }

    small_good_args = {
    'basis' : args['basis'],
    'name' : args['name'],
    'eta' : args['etam'],
    'method' : args['theorym'],
    'scratch': args['scratch']+"/cmls",
    'sign': 1,
    'og_name': file_name
    }
    small_bad_args = {
    'basis' : args['basis'],
    'name' : args['name'],
    'eta' : args['etam'],
    'method' : args['theoryw'],
    'scratch': args['scratch']+"/cmls",
    'sign': -1,
    'og_name': file_name
    }
    big_args = {
    'basis' : args['basis'],
    'name' : args['name'],
    'eta' : args['etaw'],
    'method' : args['theoryw'],
    'scratch': args['scratch']+"/cmls",
    'sign' : 1,
    'og_name': file_name
    }

    file_name = args['name'][0]
    tree = ET.parse(file_name)
    root = tree.getroot()

    dcp = kwargs.get('molecule')
    geom1py = np.asarray(dcp.geometry().to_array())
    if os.path.exists(scratchdir)==False:
        os.system('mkdir '+scratchdir)
    ofile = open(scratchdir+"/geom", "w")
    ofile.write(str(geom1py))
    ofile.close()
    basis = psi4.core.BasisSet.build(dcp, "ORBITAL", 'STO-3G')
    wfn = psi4.core.Wavefunction(dcp, basis)
    os.system("python update_cml.py "+scratchdir+"/geom "+file_name)
    frag_script.Super_Fragment(args, [big_args, small_bad_args, small_good_args])
    eg = quantum_methods.Compute_Gradient(args, root)
    true_energy = eg[1]
    true_gradient = eg[0]
    
    #dcp.update_geometry()
    geom1py = dcp.geometry().to_array()
    ofile.close()
    ofile = open(scratchdir+"/geom", "w")
    ofile.write(str(geom1py))
    grad = true_gradient
    energy = true_energy
    grad = psi4.core.Matrix.from_array(grad)
    wfn.set_gradient(grad)
    psi4.core.set_variable('CURRENT ENERGY', float(energy))
    print(true_gradient)
    print(energy)
    return (wfn)

# Integration with driver routines
psi4.driver.procedures['energy']['scf_template'] = run_scf_template
psi4.driver.procedures['gradient']['scf_template'] = run_scf_template_grad
