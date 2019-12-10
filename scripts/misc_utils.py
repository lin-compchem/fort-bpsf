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
import sys
import numpy as np


class Topology:
    """
    A class to keep track of the bonding and group topology for system
    """
    def __init__(self, numAtoms):
        self.numAtoms = numAtoms
        self.bonds = [[] for x in range(numAtoms)]
        self.atomGroups = [0] * numAtoms
        return

    def addBonds(self, a1, a2):
        """
        Add a bond to a1 and a2. This means a1 is added to a2's bonds and
        a2 is added to a1's bonds
        :param a1:
        :param a2:
        :return:
        """
        if a2 in self.bonds[a1]:
            print("Cannot add bond, atom %d is already bonded to atom %d"
                  % (a1, a2))
            raise RuntimeError
            sys.exit()
        if a1 in self.bonds[a2]:
            print("Cannot add bond, atom %d is already bonded to atom %d"
                  % (a2, a1))
            raise RuntimeError
            sys.exit()
        if a1 == a2:
            print("Cannot bond atom to itself")
            raise RuntimeError
        self.bonds[a1].append(a2)
        self.bonds[a2].append(a1)
        return

    def printBonds(self, style=None):
        pBonds = self.getUniqueBonds()
        nBonds = len(pBonds)
        bLine = "{0:9d}   {1:9d}"
        psfLine = "{0:8d}{1:8d}{2:8d}{3:8d}{4:8d}{5:8d}{6:8d}{7:8d}"
        partLine ="{0:8d}{1:8d}"
        if not style:
            print("NBonds: " + str(nBonds))
            for i in range(nBonds):
                print(bLine.format(*pBonds[i]))
        if style == 'psf':
            tBonds = []
            for i in pBonds:
                tBonds.append((i[0]+1, i[1]+1))
            pBonds = tBonds
            print("{0:8d} !NBONDS: bonds".format(nBonds))
            nLines = nBonds // 4
            b = 0
            for i in range(nLines):
                print(psfLine.format(*pBonds[b], *pBonds[b+1], *pBonds[b+2],
                                     *pBonds[b+3]))
                b += 4
            if nBonds % 4 > 0:
                lastLine = ""
                for i in range(nBonds % 4):
                    lastLine += partLine.format(*pBonds[b])
                    b += 1
                print(lastLine)
        return

    def getUniqueBonds(self):
        uniqueBonds = []
        for i in range(self.numAtoms - 1):
            for j in self.bonds[i]:
                if j > i:
                    uniqueBonds.append((i, j))
        return uniqueBonds


class PSFFile:
    def __init__(self, psfPath):
        try:
            self.ifi = open(psfPath, 'r')
        except FileNotFoundError:
            print("Cannot find psf file: " + psfPath)
            sys.exit()

    def getNumAtoms(self):
        """
        Read the number of atoms from the PSF file
        :return:
        """
        self.ifi.seek(0)
        natoms = 0
        while True:
            line = self.ifi.readline()
            if not line: break
            if "!NATOM" in line:
                try:
                    natoms = int(line.split()[0])
                except IOError:
                    print("Cannot read natoms from psf file")
                    sys.exit()
        if natoms == 0:
            print("Could not find NATOM in psf file")
            raise IOError
            sys.exit()
        return natoms

    def getBonds(self):
        numAtoms = self.getNumAtoms()
        topo = Topology(numAtoms)
        self.ifi.seek(0)
        foundNBond = False
        while True:
            line = self.ifi.readline()
            if not line: break
            if "!NBOND" in line:
                foundNBond = True
                break

        if foundNBond == False:
            print("Could not find NBonds section in psf file")
            raise IOError
            sys.exit()

        try:
            nbonds = int(line.split()[0])
        except IOError:
            print("Cannot read Nbonds")
            sys.exit()

        bondLines = nbonds // 4
        if nbonds % 4 > 0:
            bondLines += 1

        foundBonds = 0

        for l in range(bondLines):
            line = self.ifi.readline()
            words = line.split()
            nwords = len(words)
            nowBonds = nwords // 2
            for i in range(nowBonds):
                a1 = int(words[i * 2]) - 1
                a2 = int(words[i * 2 + 1]) - 1
                topo.addBonds(a1, a2)
                foundBonds += 1

        if foundBonds != nbonds:
            print("Did not find all of the bonds!")
            raise RuntimeError
            sys.exit()
        return topo


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
        '-v',
        '--verbose',
        help='Controls the level of output, use multipe v for more output',
        required=False,
        action='count',
        default=0)
    parser.add_argument(
        '-c',
        '--center',
        help='atom index to be moved to center',
        required=True,
        type=int
    )
    parser.add_argument(
        '-x',
        '--x',
        help='x length',
        required=True,
        type=float
    )
    parser.add_argument(
        '-y',
        '--y',
        help='y length',
        required=True,
        type=float
    )
    parser.add_argument(
        '-z',
        '--z',
        help='z length',
        required=True,
        type=float
    )
    parser.add_argument(
        '-d',
        '--debug',
        help='Enter debug mode',
        required=False,
        action='store_true',
        default=False)
    parser.add_argument(
        '--psf',
        help="PSF file path for reading bonds",
        required=False,
        default=None
    )
    return parser.parse_args(args)


def read_xyz(ifpath, atomic_number=False):
    """
    Reads an xyz file
    :param ifpath:
    :return: num_atoms, x, y, z as numpy int or float arrays
    """
    el_dict = {'H': 1, 'C': 6, 'N': 7, 'O': 8}
    with open(ifpath) as ifi:
        try:
            num_atoms = int(ifi.readline())
        except:
            sys.exit("Error reading n atoms")
        title = ifi.readline()
        x = np.zeros(num_atoms, dtype=float)
        y = np.zeros(num_atoms, dtype=float)
        z = np.zeros(num_atoms, dtype=float)
        atom_type = [''] * num_atoms
        for i in range(num_atoms):
            line = ifi.readline()
            words = line.split()
            atom_type[i] = words[0]
            x[i] = float(words[1])
            y[i] = float(words[2])
            z[i] = float(words[3])
            if atom_type[i] == 'OH2':
                atom_type[i] = 'O'
            elif atom_type[i] == 'H1' or atom_type[i] == 'H2'\
                    or atom_type[i] == 'H3':
                atom_type[i] = 'H'
        if atomic_number:
            for i, el in atom_type:
                atom_type[i] = el_dict[el]
            atomic_number = np.asarray(atomic_number)
    return num_atoms, x, y, z, atom_type, title


def read_xyz2(ifpath):
    """
    Reads an xyz file
    :param ifpath:
    :return: num_atoms, coords, atom_type, atom_number as numpy int or float
    :rtype: (np.ndarray(n, dtype=int16), np.ndarray(n,3, dtype=float),
    np.ndarray(n, dtype="U2"), np.ndarray(n, dtype=int8)
    """
    el_dict = {'H':1, 'C':6, 'N':7, 'O':8}
    num_dict = {1:'H', 6:'C', 7:'N', 8:'O'}
    with open(ifpath) as ifi:
        try:
            num_atoms = np.int16(ifi.readline())
        except:
            sys.exit("Error reading n atoms")
        coords = np.zeros((num_atoms,3), dtype=float)
        atom_type = np.zeros(num_atoms, dtype="U2")
        atomic_number = np.zeros(num_atoms, dtype=np.int8)
        line = ifi.readline()
        for i in range(num_atoms):
            line = ifi.readline()
            words = line.split()
            atom_type[i] = words[0]
            coords[i, 0] = float(words[1])
            coords[i, 1] = float(words[2])
            coords[i, 2] = float(words[3])
        if not atom_type[i][0].isdigit():
            for i, el in enumerate(atom_type):
                atomic_number[i] = el_dict[el]
        else:
            for i, el in enumerate(atom_type):
                atomic_number[i] = int(atom_type[i])
                atom_type[i] = num_dict[atomic_number[i]]

    return num_atoms, coords, atom_type, atomic_number


def wrap_coords(cell_x, cell_y, cell_z, x, y, z, num_atoms, center):
    """
    1. Translate all atoms by vector from atom designated by center
    to center of box defined by [[0, cell_x],[0, cell_y], [0, cell_z]]
    2. wrap all atoms outside of the above box into the box

    :param cell_x: x cell dimention float
    :param cell_y: y cell dimention float
    :param cell_z: float z cell dimension
    :param x: np float array x coords
    :param y: np float array y coords
    :param z: np float array z coords
    :param num_atoms: int
    :param center: int 0-based index of atom center
    :return: x, y, z
    """
    def wrap_dimension(x, cell_x):
        no_wraps = False
        while not no_wraps:
            no_wraps = True
            for i in range(x.size):
                if x[i] < 0:
                    x[i] += cell_x
                    no_wraps = False
                elif x[i] > cell_x:
                    x[i] -= cell_x
                    no_wraps = False
        return x



    # Find the translation vector to shift all atoms by the vector
    # from the desired centers coords to the new center
    center_a_coords = np.asarray([x[center], y[center], z[center]])
    new_center = np.asarray([cell_x/2., cell_y/2., cell_z/2.])
    translation = new_center - center_a_coords

    x[:] += translation[0]
    y[:] += translation[1]
    z[:] += translation[2]


    # Now wrap the atoms
    x = wrap_dimension(x, cell_x)
    y = wrap_dimension(y, cell_y)
    z = wrap_dimension(z, cell_z)
    return x, y, z


def print_xyz(num_atoms, x, y, z, atom_type, title=None):
    print(num_atoms)
    if title:
        print(title[:-1])
    else:
        print()
    l = "{0:2}    {1:10.6f}    {2:10.6f}    {3:10.6f}"
    for i in range(num_atoms):
        print(l.format(atom_type[i], x[i],y[i], z[i]))


def testPSFBonds():
    psf = PSFFile("/home/aduster/cp2k/h3o/proc_blyp/h3o.psf")
    bonds = psf.getBonds()
    bonds.printBonds(style='psf')

# ifpath = args.input
# cell_x = args.x
# cell_y = args.y
# cell_z = args.z
# center = args.center
# verbose = args.verbose
# debug = args.debug
#
# num_atoms, x, y, z, atom_type, title = read_xyz(ifpath)
# x, y, z = wrap_coords(cell_x, cell_y, cell_z, x, y, z, num_atoms, center)
# print_xyz(num_atoms, x, y, z, atom_type, title)

def main():
    arg_vals = None
    args = get_args(arg_vals)

if __name__ == "__main__":
    testPSFBonds()