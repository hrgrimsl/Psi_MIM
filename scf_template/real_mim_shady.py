import os
import sys
from sys import argv

script, system_name, small_eta, small_method, small_basis, big_eta, big_method, big_basis, scratchdir = argv

def Big_Bad():
    os.system("python frag_prog_gamma.py "+system_name+" "+big_eta+" "+scratchdir)
    os.system("python lastmimsy.py "+system_name+" "+big_method+" "+big_basis+" "+scratchdir)
    efile = open(scratchdir+"/final_energy","r")
    energy = float(efile.read())
    return float(energy)

def Small_Good():
    os.system("python frag_prog_gamma.py "+system_name+" "+small_eta+" "+scratchdir)
    os.system("python lastmimsy.py "+system_name+" "+small_method+" "+small_basis+" "+scratchdir)
    efile = open(scratchdir+"/final_energy","r")
    energy = float(efile.read())
    return float(energy)

def Small_Bad():
    os.system("python lastmimsy.py "+system_name+" "+big_method+" "+big_basis+" "+scratchdir)
    efile = open(scratchdir+"/final_energy","r")
    energy = float(efile.read())
    return float(energy)

big_bad = Big_Bad()
small_good = Small_Good()
small_bad = Small_Bad()

mim_energy = big_bad+small_good-small_bad
efile = open(scratchdir+"/final_energy","w")
efile.write(str(mim_energy))
efile.close()
print (mim_energy)

