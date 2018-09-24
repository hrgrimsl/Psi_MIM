import xml.etree.ElementTree as ET
import sys
from sys import argv

script, geom_file, cml_file = argv

#bohr = 0.5291772106712
bohr = 0.52917721067
tree = ET.parse(cml_file)
root = tree.getroot()
g_file = open(geom_file, 'r')
line_no = 0
for line in g_file.readlines():
    line = line.replace("[[","")
    line = line.replace("]]",'')
    line = line.replace("[","")
    line = line.replace("]",'')
    line = line.replace("\n",'')
    line = line.split(" ")
    job_done = False
    while job_done == False:
        job_done = True
        for thing in line:
            if " " in str(thing) or str(thing)=='':
                line.remove(thing)
                job_done = False
    root[0][line_no].attrib['x3']=str(float(line[0])*bohr)
    root[0][line_no].attrib['y3']=str(float(line[1])*bohr)
    root[0][line_no].attrib['z3']=str(float(line[2])*bohr)
    line_no = line_no+1
tree.write(cml_file)
