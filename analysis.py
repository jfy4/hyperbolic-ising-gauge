#!/home/jfunmuth/miniconda3/bin/python

# use taxi-cab distances between the boundary plaquettes
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsla
import sys
import os
# import matplotlib.pyplot as plt
# %matplotlib inline

def resample(data, num=1):
    """
    Computes the jackknife resampled data from a data set.
    
    Parameters
    ----------
    data : The data to be resampled.
    num  : The block size used to resample the data.
    
    Returns
    -------
    means : The resampled data of size according to 'num'.
    
    """
    if (num == 1):
        N = len(data)
        means = (np.sum(data) - data)/float(N-1)
        return means
    else:
        temp = data
        the_sum = np.sum(temp, dtype=np.float64)
        data_len = len(temp)
        k = data_len/num
        idx = [num*(i+1) for i in range(k-1)]
        blocks = np.split(temp, idx)
        sums = map(np.sum, blocks)
        means = np.array((the_sum - np.asarray(sums)) / 
                      (data_len - np.asarray(map(len, blocks))))
        return means
    
def dobin(array, n):
    data = np.asarray(array)
    size = data.shape[0]
#     data[data == 0.] = np.nan
#     split = np.array_split(data, range(n, size, n))
#     if len(split[-1]) != len(split[0]):
#         split.pop()
    split = np.array_split(data, int(size/n))
#     print(split)
    bins = np.array(list(map(np.nanmean, split)))
    return bins

def dobinspread(array, n):
    data = np.asarray(array)
    size = data.shape[0]
#     data[data == 0.] = np.nan
    split = np.array_split(data, range(n, size, n))
    bins = np.array(list(map(np.nanstd, split)))
    return bins

def unjack(vec):
    thesum = np.sum(vec)
    n = vec.shape[0]
    want = thesum - (n-1) * vec
    return want

def make_unjackknife_resamples(data, axis=0):
    """
    Makes jackknife resamples for data array.  Data
    should be of size (N_measurements x N_variables x M).
    
    Parameters
    ----------
    data       : The data to be resampled. Should be of shape
                 (N_meas x N_var x M).
    name       : The name the resamples will be stored under in the dictionary.
    block_size : The size of the blocks used in the resampling.
    
    Returns
    -------
    array : An array of the same shape as data that has been jackknifed
            resampled.

    """
    return np.apply_along_axis(unjack, axis, data)

def make_jackknife_resamples(data, axis=0):
    """
    Makes jackknife resamples for data array.  Data
    should be of size (N_measurements x N_variables x M).
    
    Parameters
    ----------
    data       : The data to be resampled. Should be of shape
                 (N_meas x N_var x M).
    name       : The name the resamples will be stored under in the dictionary.
    block_size : The size of the blocks used in the resampling.
    
    Returns
    -------
    array : An array of the same shape as data that has been jackknifed
            resampled.

    """
    return np.apply_along_axis(resample, axis, data)

def jackerror(vals):
    return np.std(vals) * np.sqrt(len(vals)-1)

def bin_scan(arr):
    binsizes = range(1,len(arr)//2)
    vals = list()
    for bs in binsizes:
        ubcrr = resample(dobin(arr, bs))
#         ubcrr = mode_resample(dobin(arr, bs))
#         print(ubcrr)
        err = jackerror(ubcrr)
        vals.append(err)
    vals = np.array(vals)
    return vals

def p_to_p(pid, neighs):
    """
    Taxi-cab distance for plaq-plaq
    """
    sphere = set([pid])
    boundary_sphere = list()
    boundary_sphere.append(sphere.copy())
    boundary = set(neighs[pid])
    while (len(boundary) != 0):
        boundary_sphere.append(boundary.copy())
        new_boundary = sum((list(neighs[point]) for point in boundary.copy()), [])
        sphere = sphere.union(boundary.copy())
        boundary = set(new_boundary).difference(sphere)
    return boundary_sphere


def domain(plaq_id, nlayers):
    """
    shelling of plaquettes up to nlayers
    """
    plaq_domain = set()
    shelling = list()
    plates = list()
    
    plaq_domain.add(plaq_id)
    bounding_edges = set(bound_n21[plaq_id])
    bounding_verts = set(bound_n20[plaq_id])
    
    shelling.append(bounding_edges.copy())
    plates.append(plaq_domain.copy())

    n = 0
    while n < nlayers:
        looping_bv = bounding_verts.copy()
        for vert in looping_bv:
            curr_plaqs = bound_n02[vert]
            # print(curr_plaqs)
            for plaq in curr_plaqs:
                if plaq in plaq_domain:
                    continue
                plaq_domain.add(plaq)
                curr_edges = set(bound_n21[plaq])
                # print(curr_edges)
                # print(bounding_edges)
                bounding_edges = bounding_edges.symmetric_difference(curr_edges)
                # print(bounding_edges)
                # exit()
                curr_verts = bound_n20[plaq]
                for v in curr_verts:
                    bounding_verts.add(v)
            bounding_verts.remove(vert)
        # print(bounding_edges)
        # exit()
        shelling.append(bounding_edges.copy())
        # print(shelling)
        # exit()
        plates.append(plaq_domain.copy())
        n += 1
    return (shelling, plates)
    



if __name__ == '__main__':

    L = "18"
    beta = np.float64(sys.argv[1])
    nlayers = 15
    max_dist = 10
    num_sources = 20
    layers = range(nlayers)
    num_configs = 10000
    
    n32 = np.genfromtxt("L" + L + "_cube_facets.txt", dtype=int)-1
    plaq_cubes = dict()
    for i in range(n32.shape[0]):
        plaqs = n32[i,:]
        for p in plaqs:
            try:
                plaq_cubes[p].append(i)
            except KeyError:
                plaq_cubes[p] = [i]
    del n32 # free up this space
        
    n21 = np.genfromtxt("L" + L + "_facet_edges.txt", dtype=int)-1
    n20 = np.genfromtxt("L" + L + "_facet_verts.txt", dtype=int)-1

    boundary_plaqs = list()
    for p, c in plaq_cubes.items():
        if len(c) == 1:
            boundary_plaqs.append(p)
        elif len(c) == 2:
            pass
        else:
            raise ValueError("too many cubes!")

    print(len(boundary_plaqs))


    bound_n21 = dict()
    bound_n12 = dict()
    plaq_neighbors = dict()
    bound_n20 = dict()
    bound_n02 = dict()
    # n21array = n21.toarray()
    # n20array = n20.toarray()
    for plaq in boundary_plaqs:
        edges = n21[plaq, :]
        bound_n21[plaq] = edges
        for edge in edges:
            try:
                bound_n12[edge].append(plaq)
            except KeyError:
                bound_n12[edge] = [plaq]
        verts = n20[plaq,:]
        bound_n20[plaq] = verts
        for vert in verts:
            try:
                bound_n02[vert].append(plaq)
            except KeyError:
                bound_n02[vert] = [plaq]
#    print(len(bound_n21), len(boundary_plaqs))
#    assert False
    for edge in bound_n12.keys():
        plaqs = bound_n12[edge]
        try:
            plaq_neighbors[plaqs[0]].append(plaqs[1])
        except KeyError:
            plaq_neighbors[plaqs[0]] = [plaqs[1]]
        try:
            plaq_neighbors[plaqs[1]].append(plaqs[0])
        except KeyError:
            plaq_neighbors[plaqs[1]] = [plaqs[0]]
        # verts = np.nonzero(n20array[plaq,:])[0]
        # bound_n20[plaq] = verts
        # for vert in verts:
        #     try:
        #         bound_n02[vert].append(plaq)
        #     except KeyError:
        #         bound_n02[vert] = [plaq]

    del n21, n20 # free up this space
    # print(plaq_neighbors)
    # rings = p_to_p(list(plaq_neighbors.keys())[0], plaq_neighbors)



    #  sort the time series
    files = os.listdir("./hb_b" + str(beta).replace('.', 'p') + "/L" + L)
    config_numbers = list()
    for file in files:
        config_numbers.append(int(file.split('_')[-1][1:-4]))
    idx = np.argsort(config_numbers)
    sorted_files = np.asarray(files)[idx]
    # print(sorted_files)
    # assert False

    num_eq_files = len(sorted_files)
    print(num_eq_files)
    cdiff = num_eq_files - num_configs

    corr_array = np.zeros((num_configs,
        num_sources, max_dist))

    wl_array = np.zeros((num_configs,
        num_sources, nlayers+1))

    dom_array = np.zeros((num_configs,
        num_sources, nlayers+1))
    
    shell_array = np.zeros((num_configs,
        num_sources, nlayers+1))

    # loop through configs
    # print(np.asarray(config_numbers)[idx][cdiff:])
    #assert False
    for ff, file in enumerate(sorted_files[cdiff:]):
        print(ff, num_configs, file)
        config_spins = np.genfromtxt("./hb_b" + str(beta).replace('.', 'p') + "/L" + L + "/" + file, dtype=int)
        # the choices of sources
        sources = np.random.choice(boundary_plaqs, size=num_sources, replace=False)
        for ii, source in enumerate(sources):
            # print(source)
            source_plaq = np.prod(config_spins[np.asarray(list(bound_n21[source]))])
            rings = p_to_p(source, plaq_neighbors)
            # loop through each distance for plaquettes
            # fill in the correlator array
            for d, ring in enumerate(rings[:max_dist]):
                # print(ring)
                # print(source_plaq, np.prod(config_spins[np.asarray(list(bound_n21[plaq]))]), d)
                corr_array[ff, ii, d] = np.mean(source_plaq*np.asarray([np.prod(config_spins[np.asarray(list(bound_n21[plaq]))])
                                                                    for plaq in ring]))
            # assert False
            # Compute the wilson loops and domains
            shells, plaq_domains = domain(source, nlayers)
            # print(plaq_domains)
            # exit()
            assert len(plaq_domains[-1]) < len(boundary_plaqs)//2
            wl_distance = list()
            domain_distance = list()
            # loop through each shell of plaquette domains
            for jj, (shell, dom) in enumerate(zip(shells, plaq_domains)):
                # print(np.asarray(list(edge_ids), dtype=int))
                wl_spins = config_spins[np.asarray(list(shell), dtype=int)]
                # print(shell)
                # assert False
                wl = np.prod(wl_spins)
                # print(wl)
                wl_array[ff, ii, jj] = wl
                dom_array[ff, ii, jj] = len(dom)
                shell_array[ff, ii, jj] = len(shell)
            # exit()
        
        # print(wl_array[0,0,:])
        # exit()

    # # print(wl_distance)
    # # print(domain_distance)
    np.save("./analL" + L + "/hb_wilson_loop_b" + str(beta).replace('.', 'p') + "_L" + L + ".npy", wl_array)
    np.save("./analL" + L + "/hb_domain_b" + str(beta).replace('.', 'p') + "_L" + L + ".npy", dom_array)
    np.save("./analL" + L + "/hb_shell_b" + str(beta).replace('.', 'p') + "_L" + L + ".npy", shell_array)
    np.save("./analL" + L + "/hb_pcorr_b" + str(beta).replace('.', 'p') + "_L" + L + ".npy", corr_array)

