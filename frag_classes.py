import numpy as np
import xml.etree.ElementTree as ET
import copy
import cProfile
import copy as cp
from frag_methods import *
from shutil import copyfile


'''
Unfragmented system
'''

class Molecule:
    def __init__(self):
        #Number of atoms
        self.n_atoms = 0
        #Atom objects
        self.atoms = []
        #Bonds of each atom
        self.bond_table = []
        #List of bond orders
        self.bond_dict = {}
        #AKA units- indivsible by whatever criteria we're using
        self.prims = []
        #Actual fragments with orders and everything
        self.frags = []
        #List of where fragments of each order can be found
        self.partition = [0]
        #List of combos we have a fragment for!
        self.combos = []
        #Dictionary of 'real' fragments for Doc's pie function
        self.file_name = ""
        self.tree = ''
        self.covrad = form_covalent_radii()


    def parse_cml(self, file_name):
        self.file_name = file_name
        tree = ET.parse(file_name)
        self.tree = tree
        root = tree.getroot()
        self.n_atoms = len(root[0])
        for atomi in range(0, self.n_atoms):
            self.atoms.append(Atom())
            self.atoms[-1].element = root[0][atomi].attrib['elementType']
            self.atoms[-1].index = atomi
            self.bond_table.append([])
            self.atoms[-1].xyz = np.asarray([ float(root[0][atomi].attrib['x3']) , float(root[0][atomi].attrib['y3']) , float(root[0][atomi].attrib['z3'])  ])
        for bondi in root[1]:
            a12 = bondi.attrib['atomRefs2'].split()
            assert(a12[0][0] == "a")
            assert(a12[1][0] == "a")
            a1 = int(a12[0][1:])-1
            a2 = int(a12[1][1:])-1
            self.bond_table[a1].append(a2)
            self.bond_table[a2].append(a1)
            self.bond_dict[(a1,a2)] = int(bondi.attrib["order"])
            self.bond_dict[(a2,a1)] = int(bondi.attrib["order"])

    def get_prims(self):
        assignment = ['no']*self.n_atoms
        for atomi in range(0, self.n_atoms):
            if assignment[atomi] == 'no':
                tmp = Primitive()
                tmp.index = len(self.prims)
                tmp.atoms = [atomi]
                self.atoms[atomi].in_prim = tmp.index
                assignment[atomi] = 'yes'
                done = False
                while done == False:
                    done = True
                    for atomj in tmp.atoms:
                        for atomk in self.bond_table[atomj]:
                            if assignment[atomk] == 'no':
                                if self.bond_dict[(atomj, atomk)]>1 or self.atoms[atomj].element == 'H' or self.atoms[atomk].element == 'H':
                                    assignment[atomk] = 'yes'
                                    done = False
                                    self.atoms[atomk].in_prim = tmp.index
                                    tmp.atoms.append(atomk)
                self.prims.append(tmp)

    def get_prim_conns(self):
        for prim in self.prims:
            conn_list = set([])
            for atomi in prim.atoms:
                conn_list = conn_list.union(set(self.bond_table[atomi]))
            prim_con_list = set([])
            for atomi in conn_list:
                if self.atoms[atomi].in_prim != prim.index:
                    prim_con_list.add(self.atoms[atomi].in_prim)
            self.prims[prim.index].conn_list = list(prim_con_list)

    def get_frags(self, eta):
         self.frags = []
         if eta!= -1:
             for prim in self.prims:
                 tmp = set([prim.index])
                 i = 0
                 append = set([])
                 while i<eta:
                     for prim2 in tmp:
                         append = append.union(self.prims[prim2].conn_list)
                     tmp = tmp.union(append)
                     i = i+1
                 self.frags.append(Fragment())
                 self.frags[-1].prims = list(tmp)

    def cull_frags(self):
        cull_list = []
        for frag1 in range(0, len(self.frags)):
            for frag2 in range(frag1+1, len(self.frags)):
                if set(self.frags[frag1].prims)==set(self.frags[frag2].prims):
                    cull_list.append(self.frags[frag2])
        for frag in cull_list:
            if frag in self.frags:
                self.frags.remove(frag)
        cull_list = []
        for frag1 in range(0, len(self.frags)):
            for frag2 in range(0, len(self.frags)):
                if self.frags[frag1]!=self.frags[frag2] and set(self.frags[frag2].prims).issubset(set(self.frags[frag1].prims)):
                    cull_list.append(self.frags[frag2])
        for frag in cull_list:
            if frag in self.frags:
                self.frags.remove(frag)
        for frag in self.frags:
            frag.atoms = []
            for prim in frag.prims:
                frag.atoms.extend(self.prims[prim].atoms)

    def close_frags(self):
        for frag in self.frags:
            Done = False
            while Done == False:
                Done = True
                neighbors = []
                for atomi in frag.atoms:
                    for atomj in self.bond_table[atomi]:
                        if atomj not in frag.atoms:
                            neighbors.append(atomj)
                restore = set([])
                for atomi in range(0, len(neighbors)):
                    for atomj in range(atomi+1, len(neighbors)):
                        if neighbors[atomi] == neighbors[atomj]:
                            restore.add(neighbors[atomi])
                            Done = False
                for atom in restore:
                    if self.atoms[atom].in_prim not in frag.prims:
                        frag.prims.append(self.atoms[atom].in_prim)
                frag.atoms = []
                for prim in frag.prims:
                    frag.atoms.extend(self.prims[prim].atoms)

    def finalize_first_frags(self):
        for fi in range(0, len(self.frags)):
            self.frags[fi].order = 1
            self.frags[fi].index = fi
            self.frags[fi].history = set([fi])

    def construct_frag_dict(self):
        self.dict = {}
        for frag in self.frags:
            frag_tuple = tuple(sorted(set(frag.prims)))
            self.dict[frag_tuple] = [frag.index, 1]
        return self.dict

    def make_frag_objects(self, final_frag_list):
        self.frags = []
        for frag in final_frag_list:
            if final_frag_list[frag][1]!=0:
                self.frags.append(Fragment())
                self.frags[-1].order = final_frag_list[frag][1]
                self.frags[-1].prims = frag
                for prim in self.frags[-1].prims:
                    self.frags[-1].atoms.extend(self.prims[prim].atoms)

    def append_meta_list(self, frag, scratchdir, index):
        order = frag.order
        dst_file = str(self.file_name).replace(".cml","")+"_"+str(index)+"_"+str(order)+".cml"
        meta = open(str(scratchdir)+"/"+str(self.file_name).replace(".cml","")+"_CML_list", "a")
        meta.write(str(dst_file)+"\n")

    def write_cml(self, frag, dst_file, args):
        copyfile(self.file_name, dst_file)
        tree = ET.parse(dst_file)
        root = tree.getroot()
        charge = 0
        mult = 1
        for target in root.findall('target'):
            id = int(target.attrib['id'])-1
            if 'X' not in root[0][id].attrib['elementType']:
                charge += int(target.attrib['charge'])
                mult += int(target.attrib['mult'])
        for atom in range(0, len(root[0])):
            if atom not in frag.atoms:
                root[0][atom].attrib['elementType'] = 'X'
            #bond = Bond(self.atoms[atom], self.atoms[atom], self)
        bonds = []
        for atom in frag.atoms:
            for atom2 in self.bond_table[atom]:
                if atom2 not in frag.atoms:
                    bonds.append(Bond(self.atoms[atom], self.atoms[atom2], self))
        for bond in bonds:
            new_H = ET.SubElement(root[0], 'atom')
            h_new = len(root[0])
            new_H.attrib['id']="a"+str(h_new)
            new_H.attrib['elementType']='H'
            new_H.attrib['x3'] = str(bond.new_xyz[0])
            new_H.attrib['y3'] = str(bond.new_xyz[1])
            new_H.attrib['z3'] = str(bond.new_xyz[2])
            new_bond = ET.SubElement(root[1], 'bond')
            new_bond.attrib['atomRefs2'] = "a"+str(bond.stay+1)+" a"+str(h_new)
            new_bond.attrib['order'] = '1'
            new_atom = ET.SubElement(root[0], 'atom')
            new_atom.attrib['elementType'] = "XR_"+str(bond.stay)+"_"+str(h_new)+"_"+str(bond.leave)+"_"+str(bond.factor)
            new_atom.attrib['id'] = 'a'+str(len(root[0]))
            new_atom.attrib['x3']=root[0][bond.leave].attrib['x3']
            new_atom.attrib['y3']=root[0][bond.leave].attrib['y3']
            new_atom.attrib['z3']=root[0][bond.leave].attrib['z3']
        if mult==1:
            child = ET.Element("theory", name=args['method'])
        else:
            child = ET.Element("theory", name='MCSCF')
        root.append(child)
        child2 = ET.Element("basis", name=args['basis'])
        root.append(child2)
        child3 = ET.Element("sign", name=str(args['sign']))
        root.append(child3)
        for i in root.findall('chargemult'):
            i.attrib['cm']=str(charge)+' '+str(mult)
        tree.write(dst_file)


    def __repr__(self):
        return "("+str(self.atoms)+")"

class Bond:
    def __init__(self, a1, a2, molecule):
        cov1 = molecule.covrad[a1.element][0]
        cov2 = molecule.covrad[a2.element][0]
        vector = a2.xyz-a1.xyz
        dist = np.linalg.norm(vector)
        h = .32
        self.stay = a1.index
        self.leave = a2.index
        self.factor = (h+cov1)/(cov1+cov2)
        self.new_xyz = self.factor*vector+a1.xyz

class Primitive:
    def __init__(self):
        self.atoms = []
        self.index = 0
        self.conn_list = []

    def __repr__(self):
        return "("+str(self.index)+")"

class Fragment:
    def __init__(self):
        self.atoms = []
        self.prims = []
        self.order = 0
        self.index = 0
        self.neighbors = {}


    def __repr__(self):
        return self.history

class Atom:
    def __init__(self):
        self.in_prim = 0
        self.element = 'X'
        self.index = 0
        self.xyz = np.asarray([])

    def __repr__(self):
        return str(self)
    def __str__(self):
        return self.element + "(" + str(self.index) + ")"
