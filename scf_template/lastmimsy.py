import os
import xml.etree.ElementTree as ET
from sys import argv
script, system_name, method, basis, scratchdir = argv
import multiprocessing

'''
Beginning of CML -> Psi4 Input
'''
def Clear_Fragment_List(psi4_fragment_list):
    frag_list = open(psi4_fragment_list, "w")
    frag_list.close()

def Get_CML_Order(line):
    newline = line.split('_')
    order = int(newline[2])
    return order

def Get_Fragment_Name(line):
    newline = line.split('_')
    fragment_name = str(newline[0])+str(newline[1])
    return fragment_name

def Update_Fragment_List(fragment_name, psi4_fragment_list):
    frag_list = open(psi4_fragment_list, "a")
    frag_list.write(str(fragment_name)+"\n")

def Create_Psi4_File(line, fragment_name, cml_order, energy_list):
    cml_file_source = str(scratchdir)+"/"+str(line)+".cml"
    cml_file_source = cml_file_source.replace("\n","")
    psi4_file = str(scratchdir)+"/"+str(fragment_name)
    tree = ET.parse(cml_file_source)
    root = tree.getroot()
    ofile = open(psi4_file, "w")
    ofile.write("set basis "+str(basis))
    ofile.write("\n")
    ofile.write("memory 10 GB")
    ofile.write("\n")
    ofile.write("molecule {")
    ofile.write("\n")
    for atom in root[0]:
        if 'X' not in atom.attrib['elementType']:
            ofile.write(" ")
            ofile.write(atom.attrib["elementType"])
            ofile.write(" ")
            ofile.write(atom.attrib["x3"])
            ofile.write(" ")
            ofile.write(atom.attrib["y3"])
            ofile.write(" ")
            ofile.write(atom.attrib["z3"])
            ofile.write("\n")
    ofile.write("}")
    ofile.write("\n")
    ofile.write("system_energy = energy(\'"+str(method)+"\')\n")
    ofile.write("order = "+str(cml_order))
    ofile.write("\n")
    ofile.write("\n")
    ofile.write("write_file = \""+str(energy_list)+"\"\n")
    ofile.write("f_out = open(write_file, \'a\') \n \n"+ "\n")
    ofile.write("f_out.write(str(order))"+"\n")
    ofile.write("f_out.write(\" \")"+"\n")
    ofile.write("f_out.write(str(system_energy))"+"\n")
    ofile.write("f_out.write(\"\\n\")"+"\n")
    ofile.write("\n")
    ofile.write("f_out.close()")
    ofile.close()


def Parse_CML_Meta_List(cml_meta_list, psi4_fragment_list, energy_list):
    with open(cml_meta_list) as f:
        Clear_Fragment_List(psi4_fragment_list)
        for line in f.readlines():
            line = line.replace(".cml","")
            cml_order = Get_CML_Order(line)
            fragment_name = Get_Fragment_Name(line)
            Update_Fragment_List(fragment_name, psi4_fragment_list)
            Create_Psi4_File(line, fragment_name, cml_order, energy_list)
'''
End of CML -> Psi4 Input
'''


'''
Beginning of Psi4-Calling
'''
def Clear_Energy_List(energy_list):
    clearfile = str(energy_list)
    clearfile2 = open(clearfile, "w")
    clearfile2.close()
def Run_Job(argline):
    os.system(argline)
def Make_Energy_List(fragment_list, energy_list):
    Clear_Energy_List(energy_list)
    with open(fragment_list) as f:
        jobs = []
        for line in f.readlines():
            line_literal = str(scratchdir)+"/"+str(line)
            argline = 'psi4 '+ line_literal
            if __name__ == '__main__':
                process = multiprocessing.Process(target = Run_Job, args = (argline,))
                jobs.append(process)
        for j in jobs:
            j.start()
            j.join()
'''
End of Psi4-Calling
'''



'''
Beginning of Energy-Summing
'''
def Get_Order(line):
    newline = line.split()
    return int(newline[0])

def Get_Energy(line):
    newline = line.split()
    return float(newline[1])

def Combine_Energies(energy_list):
    energy_sum = 0
    with open(energy_list) as f:
        for line in f.readlines():
            order = Get_Order(line)
            energy = Get_Energy(line)
            energy_sum = energy_sum + energy*(-1)**(order+1)
    return energy_sum

'''
End of Energy-Summing
'''


#Procedural List
cml_list = str(scratchdir)+"/"+str(system_name)+"_CML_list"
psi4_frag_list = str(scratchdir)+"/"+str(system_name)+"_fragment_list"
energy_list = str(scratchdir)+"/"+str(system_name)+"_energies"
Parse_CML_Meta_List(cml_list, psi4_frag_list, energy_list)
Make_Energy_List(psi4_frag_list, energy_list)
eout = open(scratchdir+"/final_energy","w")
eout.write(str(Combine_Energies(energy_list)))

