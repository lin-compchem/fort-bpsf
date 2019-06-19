#!/usr/bin/env python
"""
Create a fort-bpsf file using human-readable input
"""
__author__ = 'Adam Duster'
__copyright__ = ''
__credits__ = ['Adam Duster']
__license__ = 'CC-BY-SA'
__version__ = '0.1'
__email__ = 'adam.duster@ucdenver.edu'
__status__ = 'Development'

import argparse
import numpy as np
import sys


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
        required=False)
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
    parser.add_argument(
        '-e',
        '--example',
        help='Print example input file, then exit. This file can be used with'
             ' -i to create an input file for the fortran program',
        default=None,
        action='store_true'
    )
    return parser.parse_args(args)




class Element:
    """
    Parameters for each element
    """
    def __init__(self, atm_num):
        self.atm_num = int(atm_num)
    num_bonds = 0
    bond_type = []
    eta = []
    rs = []
    bond_len = []
    max_bond = 0

    num_angles = 0
    ang_len = []
    zeta = []
    lam = []
    e2 = []
    max_angle = 0

    def add_bond(self, bond_type, etas, rs):
        """
        Add the list of etas and rs to the structure
        :param bond_type:  int
        :param etas: list of floats
        :param rs: list of floats
        :return:
        """
        num_bonds += 1
        self.bond_type.append(bond_type)
        self.etas.append(etas)
        self.rs.append(rs)

        if len(etas) != len(rs):
            print("Error eta and rs have different length")
            raise
        if (not etas) or (not rs):
            print("Error initializing etas rs")
            raise
        self.bond_len.append(len(etas))

    def add_angle(self, ang_type, zeta, lam, e2):
        """
        Add the angle information
        :param ang_type: int
        :param zeta: list of floats
        :param lam: list of floats
        :param e2: list of floats
        :return:
        """
        num_angles += 1
        self.ang_type.append(ang_type)
        self.zeta.append(zeta)
        self.lam.append(lam)
        self.e2.append(e2)
        self.ang_leg.append(len(zeta))

    def count_data(self):
        self.max_bond = max(bond_len)
        self.max_angle = max(ang_len)


def read_element(ifi, el_num):
    """
    Read data into the element class
    :param ifi:
    :param el_num: integer atomic number
    :return:
    """
    e = Element(el_num)
    while True:
        l = rline(ifi)
        if l[:4] == 'bond':
            w = l.split()
            bt = int(w[1])
            eta, rs = rlist(ifi)
            e.add_bond(bt, eta, rs)
        elif l[:5] == 'angle':
            w = l.split()
            at = int(w[1])
            zeta, lam, e2 = rlist(ifi)
            e.add_angle(at, zeta, lam, e2)
        elif l[:3] == 'end':
            break
    e.count_data()
    return e


def rline(ifi):
    """
    Read a file skipping blank and comment lines
    return the first found line lowercase
    :param ifi:
    :return:
    """
    while True:
        line = ifi.readline()
        if not line:
            raise
        l = line.lower()
        l = l.lstrip()
        if not l:
            continue
        if l[0] == '#':
            continue
        return l


def rlist(ifi):
    """
    Read a file, scraping data from the columns into lists.
    Stop when end is reached
    return the lists
    :param ifi: File
    :return: list of lists
    """
    line = rline(ifi)
    num_l = len(line)
    if num_l <= 0:
        print("Error reading list")
        raise
    w = line.split()
    dat = []
    for i in num_l:
        dat.append([])
        dat[i].append(float(w[i]))
    while True:
        line = rline(ifi)
        if 'end' in line[:3]:
            break
        w = line.split()
        for i in num_l:
            dat[i].append(float(w[i]))
    return dat


def find_max_bas(els, pars):
    """
    Find the maximum number of bond and angle functions in the els list.
    Add this number
    :param els:
    :param pars:
    :return:
    """
    max_angle = 0
    max_bond = 0
    for e in els:
        if e.max_angle > max_angle:
            max_angle = e.max_angle
        if e.max_bond > max_bond:
            max_bond = e.max_bond
    pars['max_angle'] = max_angle
    pars['max_bond'] = max_bond


def default_pars():
    """
    Function which returns the default parameters for the input file.
    -------
    pars: dict
        The default dictionary of parameters
    """
    pars = {
        'max_angle': 150,
        'max_bonds': 150,
        'calc_basis': True,
        'rc': None
    }
    return pars


def print_example():
    """
    Print an example input file to use with this script.

    It is the recommended BPSF from the paper

    Returns
    -------

    TODO
    ----
    Find citation for paper

    """

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

    pline = '{:} {:}'
    elline = 'element {:2d}'
    bdline = 'bond {:2d}'
    bb = '    {:10f}     {:10f}'
    aa = '    {:10f}     {:10f}    {:10f}'
    agline = 'angle {:2d} {:2d}'
    print('rc ', Rc)
    for i in range(num_els):
        print(elline.format(els[i]))
        for j, t in enumerate(rad_types[i]):
            print(bdline, rad_types[j])
            for b in range(rad_rs[i][j].size):
                print(bb.format(rad_rs[i][j][b], rad_etas[i][j][b]))
            print('end bond')
        for j, t in enumerate(ang_types[i]):
            print(agline.format(*ang_types[i][j]))
            for a in range(zles[i,j,:].size):
                print(aa.format(*zles[i,j,:]))
            print('end angle')
        print('end element')
    sys.exit()


def main():

    ## Vars
    arg_vals = None
    args = get_args(arg_vals)
    if args.example:
        print_example()

    pars = default_pars()

    els = []

    ifi = open(args.input, 'r')
    while True:
        line = ifi.readline()
        if not line: break
        l = line.lstrip()
        if not l:
            continue
        if l[0] == '#':
            continue

        w = l.lower().split()
        kw = w[0]
        if kw == 'rc':
            try:
                pars[kw] = float(w[0])
                pars['calc_rc'] = False
            except:
                print("Error setting Rc")
                raise
        elif kw == 'element':
            num_els += 1
            try:
                el_num = int(w[1])
            except:
                print("Error reading element number after element keyword")
                raise
            els.append(read_element(ifi, el_num))
        elif kw == 'max_basis':
            try:
                pars[kw] = int(w[0])
                pars['calc_basis'] = False
            except:
                print("Error setting Rc")
                raise
    ifi.close()
    if pars['calc_basis']:
        find_max_bas(els, pars)

    print_ifi(pars, els)

if __name__ == '__main__':
    main()