import random
import numpy as np
import copy as cp
import cython_pie

# cython: profile=True

def add_up_atoms(frag_list, n_atoms):
    mult = [0]*n_atoms
    for frag in frag_list:
        for atom in frag:
            mult[atom] += frag_list[frag][1]
    print(" Number of times each atom is counted:")
    print(mult)




def get_next_layer(next_layer, final_frag_list, sign):
    sign2 = -1*sign
    for f2 in next_layer:
        next_layer2 = {}
        for f3 in next_layer:
            if next_layer[f2][0] >= next_layer[f3][0]:
                continue
            f23 = tuple(sorted(set(f2).intersection(f3)))
            if len(f23) > 0:
                if f23 in next_layer2:
                    next_layer2[f23][1] += 1
                else:
                    next_layer2[f23] = [len(next_layer2), 1]
        for f2 in next_layer2:
            if next_layer2[f2][1] == 0:
                continue
            if f2 in final_frag_list:
                final_frag_list[f2][1] += 1*sign2
            else:
                final_frag_list[f2] = [len(final_frag_list), 1*sign2]
        get_next_layer(next_layer2, final_frag_list, sign2)


def pie(primary_frags):
    """
    This function expects a dictionary of the following form:
        [(3,6,7,8)] = [frag#, frag coeff]
             ^           ^        ^
             |           |        |
             a           b        c
        a = sorted tuple of primitives. i.e,  tuple(sorted(set(f)))
        b = index of fragment
        c = fragment coefficient. I assume these will all start out == 1

    """

    final_frag_list = cp.deepcopy(primary_frags)
    sign = -1
    for f1 in primary_frags:
        next_layer = {}
        layer = 0
        frags_curr = primary_frags
        for f2 in frags_curr:
            if frags_curr[f1][0] >= frags_curr[f2][0]:
                continue
            f12 = tuple(sorted(set(f1).intersection(f2)))
            if len(f12) > 0:
                if f12 in next_layer:
                    next_layer[f12][1] += 1
                else:
                    next_layer[f12] = [len(next_layer), 1]

        for f2 in next_layer:
            if f2 in final_frag_list:
                final_frag_list[f2][1] += 1*sign
            else:
                final_frag_list[f2] = [len(final_frag_list), 1*sign]

        get_next_layer(next_layer, final_frag_list, sign)
        #cython_pie.cython_get_next_layer(next_layer, final_frag_list, sign)

    return final_frag_list

#########################################################################
'''
n_atoms = 1000
random.seed('2')
mult = []
inv_frags = []
for i in range(n_atoms):
    mult.append(0)
    inv_frags.append([])

frag_size = 35
primary_frags = {}
for a in range(n_atoms-frag_size+1):
    frag = tuple(sorted(random.sample(range(n_atoms), frag_size)))
    #frag = tuple(range(a,a+frag_size))
    if frag in primary_frags:
        primary_frags[frag][1] = 1
    else:
        primary_frags[frag] = [len(primary_frags), 1]

    for ai in frag:
        mult[ai] += 1
        inv_frags[ai].append(a)

for f in primary_frags:
    print(f, primary_frags[f])

print("mult:", mult)
print(inv_frags)


final_frag_list = cython_pie.cython_pie(primary_frags)
#final_frag_list = pie(primary_frags)
'''
