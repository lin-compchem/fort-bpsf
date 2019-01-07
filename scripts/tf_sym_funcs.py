#!/usr/bin/env python
""" 
"""
__author__ = 'Adam Duster'
__copyright__ = ''
__credits__ = ['Adam Duster']
__license__ = 'CC-BY-SA'
__version__ = '0.1'
__email__ = 'adam.duster@ucdenver.edu'
__status__ = 'Development'

import argparse
import h5py as h5
import numpy as np
import math
import itertools
import time

class BasisSettings():
    def __init__(self, Rc, Rs, etas, zetas, lams, eprimes):
        pass
    pass

def get_args(args=None):
    """ This is written as a default funtion to put at beginning of all Python
    scripts which require command line arguments. This uses the argparse module
    which must be declared in the main program to ensure that the object is able
    to be used by the caller
    --Adam Duster 21 June 2017
    """
    parser = argparse.ArgumentParser(description='see header of python script')
    parser.add_argument(
        '-i',
        '--input',
        help='Input file name',
        required=True)
    parser.add_argument(
        '-o',
        '--output_path',
        help="The folder for the h5py files",
        required=True,
        default="./"
    )
    parser.add_argument(
        '--out_prefix',
        help="Prefix to add for h5py files",
        required=False,
        default=''
    )
    parser.add_argument(
        '-v',
        '--verbose',
        help='Controls the level of output, use multipe v for more output',
        required=False,
        action='count',
        default=0)
    parser.add_argument(
        '-d',
        '--debug',
        help='Enter debug mode',
        required=False,
        action='store_true',
        default=False)
    return parser.parse_args(args)


def getG2RadElement(fc, distance, eta, R_s, debug=False):
    """
    Radial distance vector G2 as described by Behler,
    J. Chem. Phys. 134, 074106 (2011)
    """
    # if debug:
    #     a = np.exp(-eta * np.square(distance - R_s)) * fc
    #     print("eta", eta, "Rs", R_s)
    #     print("fcs", fc)
    #     print("dists", distance)
    #     print(a)
    #     print(np.sum(a))
    #raise
    return np.sum(np.exp(-eta * np.square(distance - R_s)) * fc)


def getG4AngBasis(thetas, fcs, summed_dists, zeta, lam, eta, debug=False):
    """
    Angular symmetry scalar G4 as described by Behler
    J. Chem. Phys. 134, 074106 (2011)
    """
    # fs = np.exp(-1 * eta * summed_dists) ** 2
    # fs *= 1
    # if debug:
    #     if (eta == 0.001) and (zeta == 1.0) and (lam == -1.0):
    #         radfilter = np.exp(-eta * (summed_dists) ** 2)
    #         gangijk = 2 ** (1.0 - zeta) * (
    #                     1.0 + lam * thetas) ** zeta * radfilter
    #         print(zeta, lam, eta, 'zeta, lambda, eta')
    #         print(thetas, 'cosangles')
    #         print(fcs, 'fcs')
    #         print(summed_dists, 'sum_dists')
    #         print(radfilter, 'radfilters')
    #         print(gangijk, 'gangijk')
    #         gangijk *= fcs
    #         print(gangijk, 'after fc')
    #         gangijk = np.sum(gangijk)
    #         print(np.sum(gangijk), 'final')
    #         if not np.isclose(gangijk, 2 ** (1.0 - zeta) *
    #                                    np.sum((1 + lam * thetas) ** zeta *
    #                                           np.exp(-1 * eta * (
    #                                           summed_dists) ** 2)
    #                                           * fcs)):
    #             print('function output:')
    #             print(2 ** (1.0 - zeta) * np.sum((1 + lam * thetas) ** zeta *
    #                                              np.exp(-1 * eta * (
    #                                                  summed_dists) ** 2)
    #                                              * fcs))
    #             raise

    return 2 ** (1.0 - zeta) * np.sum((1 + lam * thetas) ** zeta *
                                      np.exp(-1 * eta * summed_dists**2) * fcs)


def applyCutoff(pair_dists, cutoff=8.):
    """
    the cutoff vector fc(rij) as described by

    dist_vec values and cutoff_vec values are validated
    """
    cutoff_vec = 0.5 * (np.cos(math.pi * pair_dists / cutoff) + 1)
    cutoff_vec[pair_dists > cutoff] = 0.
    return cutoff_vec


def main():
    test = True
    arg_vals = None
    if test:
        arg_vals = ['-i', '../test_files/2_water_clusters.h5',
                '-o', '../test_files/reference.h5',
                    '-v']
#     if test:
#         arg_vals = ['-i', '/ml/mm_md/partitions/h5s/mm_md-cartgeom.h5',
#                 '-o', '/ml/mm_md/partitions/h5s/mm_md-symfuncs.h5',
#                     '-v']
    args = get_args(arg_vals)

#############   Set PARAMETERS   ###############
    Rc = 8.0
    Rs = np.linspace(0.8, Rc, 24)
    etas = 0.5 * (5. / Rs) ** 2
    zetas = np.array([1., 4., 16])
    lam = np.array([-1, 1])
    eprime = np.array([0.001, 0.01, 0.05])
    els = np.asarray([1,8])
    num_els = els.size
    atm_names = ('H', 'O')
    rad_types = ((1, 8), (1, 8))
    rad_etas = ((etas, etas), (etas, etas[2:]))
    rad_rs = ((Rs, Rs), (Rs, Rs[2:]))
    rad_length = (rad_rs[0][0].size + rad_rs[0][1].size,
                  rad_rs[1][0].size + rad_rs[1][1].size)
    ang_types = (((1, 1), (1, 8)), ((8, 8), (1, 8), (1, 1)))
    ang_length = (len(ang_types[0]) * zetas.size * lam.size * eprime.size,
                  len(ang_types[1]) * zetas.size * lam.size * eprime.size)
    sel_dict = {
        'H': 0,
        'O': 1,
        1: 0,
        8: 1
    }
    zles = np.zeros((num_els, zetas.size * lam.size * eprime.size, 3),
                    dtype=float)  # zeta lambda eprimes
    # TODO: change name now that numbering has changed
    for j in range(num_els):
        i = 0
        for e in eprime:
            for z in zetas:
                for l in lam:
                    zles[j, i, 0] = z
                    zles[j, i, 1] = l
                    zles[j, i, 2] = e
                    i += 1
###############################################


    ifi = h5.File(args.input, "r")
    ofi = h5.File(args.output_path, 'w')
    "Read the information from the input file"
    print(ifi["atomic_number"])

    # Find the dimensions of the output matrix
    els = np.unique(ifi["atomic_number"][:])
    els = els[els != 0]
    num_els = els.size

    num_h_sym_vecs = np.count_nonzero(ifi["atomic_number"][:] == 1)
    num_o_sym_vecs = np.count_nonzero(ifi["atomic_number"][:] == 8)
    num_mols = ifi["num_atoms"].size

    # Find the length of the sym vec
    if args.verbose:
        print("After parsing input hdf5 file")
        print("Found %d elements" % num_els)
        print("Elements:", *els)
        print("Parsing %d molecules" % num_mols)
        print("Num sym_vecs for each element: %d, %d" % (num_h_sym_vecs, num_o_sym_vecs))

    chunk_size = 10000
    chunk_buffer = chunk_size // 100
    print("Creating ofi disk h5py stuf")
    coords = np.empty(ifi['cartesian_coords'].shape, dtype=float)
    h_rad = np.empty((num_h_sym_vecs, rad_length[0]), dtype=float)
    o_rad = np.empty((num_o_sym_vecs, rad_length[1]), dtype=float)
    h_ang = np.empty((num_h_sym_vecs, ang_length[0]), dtype=float)
    o_ang = np.empty((num_o_sym_vecs, ang_length[1]), dtype=float)
    h_mol = np.empty(num_h_sym_vecs, dtype=int)
    o_mol = np.empty(num_o_sym_vecs, dtype=int)
    num_atoms = np.empty(num_mols, dtype=int)
    atom_nums = np.empty((ifi["atomic_number"].shape), dtype=ifi["atomic_number"].dtype)
    h_inds = np.empty((num_mols, 2), dtype=int)
    o_inds = np.empty((num_mols, 2), dtype=int)
    #h_ang_inds = np.empty((num_mols, 2), dtype=int)
    #o_ang_inds = np.empty((num_mols, 2), dtype=int)
    inds = [h_inds, o_inds]

    disk_h_rad = ofi.create_dataset("h_radial_sym_funcs", shape=(num_h_sym_vecs, rad_length[0]), dtype=float)
    disk_o_rad = ofi.create_dataset("o_radial_sym_funcs", shape=(num_o_sym_vecs, rad_length[1]), dtype=float)
    disk_h_ang = ofi.create_dataset("h_angular_sym_funcs",shape=(num_h_sym_vecs, ang_length[0]), dtype=float)
    disk_o_ang = ofi.create_dataset("o_angular_sym_funcs",shape=(num_o_sym_vecs, ang_length[1]), dtype=float)
    #disk_h_mol = ofi.create_dataset("h_containing_moleucle",shape=(num_h_sym_vecs,), dtype=int)
    #disk_o_mol = ofi.create_dataset("o_containing_moleucle",shape=(num_o_sym_vecs,), dtype=int)
    disk_num_atoms =  ofi.create_dataset("atoms_per_molecule", shape=(num_mols,), dtype=int)
    disk_h_inds = ofi.create_dataset("h_mol_index", shape=(num_mols, 2), dtype=int)
    disk_o_inds = ofi.create_dataset("o_mol_index", shape=(num_mols, 2), dtype=int)
    rad_bas = [h_rad, o_rad]
    ang_bas = [h_ang, o_ang]
    mol = 0
    ss = np.zeros(num_els, dtype=int)
    while mol < num_mols:
        print("Beginning symfunc creation!")
        start = mol
        if num_mols - mol < chunk_size:
            chunk_size = num_mols - mol
        end = chunk_size + start
        print("reading natoms and atomic numbers")
        num_atoms[:] = ifi['num_atoms'][:]
        atom_nums[:] = ifi['atomic_number'][:]

        coords = ifi['cartesian_coords'][:]
        stime = time.time()
        print("beginning major loop")
        for chunk in range(num_mols):
            if(chunk % 100 == 0):
                print("at mol: {0} ;time for 100 mols : {1}".format(chunk,
                          stime - time.time()) )
                stime = time.time()
            natoms = num_atoms[mol]
            sel = []
            sel.append(np.argwhere(atom_nums[mol,:natoms] == 1))
            sel.append(np.argwhere(atom_nums[mol,:natoms] == 8))
            for el in range(num_els):
                if mol == 0:
                    inds[el][mol, 0] = 0
                else:
                    inds[el][mol, 0] = inds[el][mol - 1, 1]
                inds[el][mol, 1] = inds[el][mol, 0] + sel[el].size
                ss[el] = inds[el][mol, 0]
            mat_dim = num_atoms[mol]
            up_inds = np.triu_indices(mat_dim, 1)
            lo_inds = np.tril_indices(mat_dim, -1)
            nvecs = len(up_inds[0])
            vecs = np.empty((nvecs, 3), dtype=float)
            vec_mat = np.zeros((natoms, natoms, 3), dtype=float)
            k = 0

            for i in range(natoms - 1):
                j = i+1
                m = natoms - j
                l = m + k
                vecs[k:l] = np.subtract(coords[chunk, j:natoms,:], coords[chunk, i, :])
                k += m
            vec_mat[up_inds] = vecs
            for i in range(natoms - 1):
                for j in range(i + 1, natoms):
                    vec_mat[j, i] = -1 * vec_mat[i, j]

            pdists = np.linalg.norm(vecs, axis=1)
            cuts = applyCutoff(pdists)
            pair_dists = np.zeros((mat_dim, mat_dim), dtype=float)
            cutoffs = np.zeros((mat_dim, mat_dim), dtype=float)
            pair_dists[up_inds] = pdists
            pair_dists[lo_inds] = pair_dists.T[lo_inds]
            cutoffs[up_inds] = cuts
            cutoffs[lo_inds] = cutoffs.T[lo_inds]


            # # DON"T FORGET TO CHANGE THIS!!! TODO: CHANGE THIS
            # for i in range(natoms):
            #     for j in range(natoms):
            #         if i == j:
            #             continue
            #         cutoffs[i, j] = 1.0

            for j in range(num_els):
                for m, me in enumerate(sel[j]):
                    my_el = 0
                    for r, rtype in enumerate(rad_types[j]):
                        you = sel[sel_dict[rtype]]
                        for e in range(rad_etas[j][r].size):
                            try:
                                rad_bas[j][ss[j] + m, my_el] = getG2RadElement(
                                cutoffs[me, you], pair_dists[me, you],
                                rad_etas[j][r][e], rad_rs[j][r][e], False)
                            except:
                                print("Fatal Error rad")
                                raise
                            my_el += 1
####################### ANGS
            b = time.time()
            ang_theta = np.zeros((natoms, natoms, natoms), dtype=np.float)
            ang_sum_dists = np.zeros((natoms, natoms, natoms), dtype=np.float)
            ang_fc = np.zeros((natoms, natoms, natoms), dtype=np.float)
            #print("initial setup:", time.time() - b)
            b = time.time()
            for i in range(natoms):
                for j in range(0, natoms - 1):
                    ang_theta[i, j, j + 1:] = np.matmul(vec_mat[i, j],
                                                        vec_mat[i, j + 1:].T)
            for i in range(natoms):
                for j in range(0, natoms - 1):
                    for k in range(j + 1, natoms):
                        # try:
                            ang_theta[i, j, k] /= pair_dists[i, j] * pair_dists[
                                i, k]
                            ang_sum_dists[i, j, k] = pair_dists[i, j] + \
                                                     pair_dists[i, k] + \
                                                     pair_dists[j, k]
                            ang_fc[i, j, k] = cutoffs[i, j] * cutoffs[i, k] * \
                                              cutoffs[j, k]
                        #                 if (i == 0) and (j == 3) and (k == 6):
                        #                     print('ij 0-3', pair_dists[i,j])
                        #                     print('ij 0-6', pair_dists[i,k])
                        #                     print('ij 6-3', pair_dists[j,k])
                        #                     print(ang_sum_dists[i,j,k])
                        #                     raise
                        # except:
                        #     if pair_dists[i, j] == 0 or pair_dists[i, k] == 0:
                        #         ang_theta[i, j, k] = 0
                        #         ang_sum_dists[i, j, k] = 0
                        #         ang_fc[i, j, k] = 0
                        #     else:
                        #         print("unspecified error")
                        #         raise
                ang_theta[i, i, :] = 0.
                ang_theta[i, :, i] = 0.
                np.fill_diagonal(ang_theta, 0.)
                ang_fc[i, i, :] = 0
                ang_fc[i, :, i] = 0
                np.fill_diagonal(ang_fc, 0.)

                ang_theta[i][lo_inds] = ang_theta[i].T[lo_inds]
                ang_sum_dists[i][lo_inds] = ang_sum_dists[i].T[lo_inds]
                ang_fc[i][lo_inds] = ang_fc[i].T[lo_inds]
                # pass
            #                 print(i, j, k, "i, j, k")
            #                 print(vec_mat[i,j])
            #                 print(vec_mat[i,k])
            #                 print(pair_dists[i,j])
            #                 print(pair_dists[i,k])
            # raise
            #print("setup time:", b - time.time())
            b = time.time()

            for i in range(num_els):
                aa = 0
                for a, atype in enumerate(ang_types[i]):
                    #print(a, atype)
                    #print(a, atype[0])

                    if atype[0] != atype[1]:
                        xind, yind = np.meshgrid(sel[sel_dict[atype[0]]],
                                                 sel[sel_dict[atype[1]]])
                    else:
                        s = sel[sel_dict[atype[0]]]
                        combs = np.asarray(
                            list(itertools.combinations(s, 2))).reshape(-1, 2)
                        xind, yind = combs[:, 0], combs[:, 1]
                    for p, pars in enumerate(zles[i]):
                        for m, me in enumerate(sel[i]):
                            # print()
                            try:
                                #                     if (atype[1] == 1) and (me == 0) and (atype[0] == 1):
                                #                         debug = True
                                #                         print('here')
                                #                         print(me)
                                #                         for lll, x in enumerate(xind):
                                #                             print(x, yind[lll])
                                #                             print('hi')
                                #                         raise
                                #                     else:
                                #                         debug = False
                                debug = False
                                ang_bas[i][ss[i] + m, aa] = getG4AngBasis(
                                    ang_theta[me, xind, yind],
                                    ang_fc[me, xind, yind],
                                    ang_sum_dists[me, xind, yind], pars[0],
                                    pars[1], pars[2], debug=debug)
                            except:
                                print('fatal error calculating abas')
                                raise
                        #                 if m == 1 and pars[0] == 1.0 and pars[2] == 0.001 and pars[1] == -1: raise
                        aa += 1
            mol += 1
    disk_h_rad[:] = h_rad[:]
    disk_o_rad[:] = o_rad[:]
    disk_h_ang[:] = h_ang[:]
    disk_o_ang[:] = o_ang[:]
    #disk_h_mol[:] = h_mol[:]
    #disk_o_mol[:] = o_mol[:]
    disk_num_atoms[:] = ifi["num_atoms"][:]
    disk_h_inds[:] = h_inds[:]
    disk_o_inds[:] = o_inds[:]

    ifi.close()
    ofi.close()
    #print("Check")
    #print("Check")
    return

if __name__ == "__main__":
    main()
