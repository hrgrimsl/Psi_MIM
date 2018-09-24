cimport cpython
cimport cpython.dict
#cython: boundscheck=False
#cython: nonecheck=True


cdef cython_get_next_layer(dict next_layer, dict final_frag_list, int sign):
    cdef int sign2 = -1*sign
    cdef dict next_layer2 = {}
    for f2 in next_layer:
        next_layer2 = {}
        for f3 in next_layer:
            if next_layer[f2][0] >= next_layer[f3][0]:
                continue
            f23 = set(f2).intersection(f3)
            if len(f23) > 0:
                f23 = tuple(sorted(f23))
                if f23 in next_layer2:
                    next_layer2[f23][1] += 1
                else:
                    next_layer2[f23] = [len(next_layer2), 1]
        for f2 in next_layer2:
            if next_layer2[f2][1] == 0:
                continue
            if f2 in final_frag_list:
                final_frag_list[f2][1] += 1*sign2
                #if final_frag_list[f2][1] == 0:
                #    del final_frag_list[f2]
            else:
                final_frag_list[f2] = [len(final_frag_list), 1*sign2]
        cython_get_next_layer(next_layer2, final_frag_list, sign2)


cpdef cython_pie(primary_frags):
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

    cdef dict final_frag_list = {}
    cdef list tmp2 = []
    cdef tuple tmp1 = ()
    for f1 in primary_frags:
        tmp1 = f1
        tmp2 = primary_frags[f1]
        final_frag_list[tmp1] = tmp2

    cdef sign = -1
    cdef dict next_layer = {}
    for f1 in primary_frags:
        next_layer = {}
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

        cython_get_next_layer(next_layer, final_frag_list, sign)

    return final_frag_list
