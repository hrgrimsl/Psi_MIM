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
np.set_printoptions(precision=11)

os.system('cd ..')
os.system('pwd')
si_file = './SI.txt'
sifile = open(str(si_file), 'r')
lines = sifile.readlines()
for i in range(0, len(lines)):
    lines[i] = (str(lines[i])).strip('\n')
print(lines)
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
if os.path.exists(scratchdir)==False:
    os.system('mkdir '+scratchdir)


def run_scf_template(name, **kwargs):
    '''
    #dcp = kwargs.get('molecule')
    #geom1py = dcp.geometry().to_array()
    #ofile = open("geom", "w")
    #ofile.write(str(geom1py))
    os.system("python update_cml.py geom "+file_name)
    os.system('python energy_wrapper.py '+file_name+" -nm "+lil_eta+" -tm "+lil_met+" -b "+lil_bas+" -nw "+big_eta+" -tw "+big_met+" -s "+scratchdir)
    efile = open(str(scratchdir)+"/amen","r")
    energy = float(efile.read())
    dcp = kwargs.get('molecule')
    #dcp.update_geometry()
    geom1py = dcp.geometry().to_array()
    ofile = open("geom", "w")
    ofile.write(str(geom1py))
    os.system("python update_cml.py geom "+file_name)
    psi4.core.set_variable('CURRENT ENERGY', float(energy))
    return energy
    '''
    return psi4.core.get_variable('CURRENT ENERGY')

def run_scf_template_grad(name, **kwargs):
    dcp = kwargs.get('molecule')
    geom1py = dcp.geometry().to_array()
    ofile = open(scratchdir+"/geom", "w")
    ofile.write(str(geom1py))
    ofile.close()
    #ofile = open(scratchdir+"/geom","r")
    # new skeleton wavefunction w/mol, highest-SCF basis (just to choose one), & not energy
    basis = psi4.core.BasisSet.build(dcp, "ORBITAL", 'STO-3G')
    wfn = psi4.core.Wavefunction(dcp, basis)
    os.system("cd /home/harper/PIWS")
    os.system("python update_cml.py "+scratchdir+"/geom "+file_name)
    os.system('python grad_wrapper.py '+system_name+".cml -nm "+lil_eta+" -tm "+lil_met+" -b "+lil_bas+" -nw "+big_eta+" -tw "+big_met+" -s "+scratchdir)
    gfile = open(scratchdir+"/timshel","r")
    grad = []
    for line in gfile.readlines():
        line = line.replace("\n",'')
        line = line.replace('[','')
        line = line.replace(']','')
        line = line.split()
        grad.append([float(line[0]), float(line[1]), float(line[2])])
    grad = np.asarray(grad, float)
    #dcp.update_geometry()
    geom1py = dcp.geometry().to_array()
    ofile.close()
    ofile = open(scratchdir+"/geom", "w")
    ofile.write(str(geom1py))
    os.system("python update_cml.py "+scratchdir+"/geom "+file_name)
    grad = psi4.core.Matrix.from_array(grad)
    wfn.set_gradient(grad)
    gfile.close()

    efile = open(str(scratchdir)+"/amen","r")
    energy = float(efile.read())
    psi4.core.set_variable('CURRENT ENERGY', float(energy))
    #run_scf_template(name, **kwargs)
    return (wfn)

# Integration with driver routines
psi4.driver.procedures['energy']['scf_template'] = run_scf_template
psi4.driver.procedures['gradient']['scf_template'] = run_scf_template_grad
