#!/usr/bin/env python3

# Usage:
# xyz2lmpbond.py <atom_type> <mass> <molecule_type> <starting_atom_ID> < data.xyz > data.out
# Convert a XYZ file with a single atom type to a LAMMPS data file using the 'bond' atom style.

import sys
import argparse
from enum import Enum

# A 'bond' atom style atom
class AtomBond:
    # initialize from xyz file line and some known values
    def __init__(self, xyz_line, atom_id, mol_id, atom_type):
        contents = xyz_line.split()
        num_values = len(contents)
        if num_values != 4:
            raise Exception("A xyz line is expected to have 4 values: {0:d} found." \
            "\nLine: {1:s}".format(num_values, xyz_line))
        self.atom_id = atom_id
        self.mol_id = mol_id
        self.atom_type = atom_type
        self.x = float(contents[1])
        self.y = float(contents[2])
        self.z = float(contents[3])
        self.nx = 0
        self.ny = 0
        self.nz = 0

    def __str__(self):
        return "{0:10d} {1:7d} {2:7d} {3:9.4f} {4:9.4f} {5:9.4f} {6:3d} {7:3d} {8:3d}".format(self.atom_id,
                                                                              self.mol_id,
                                                                              self.atom_type,
                                                                              self.x,
                                                                              self.y,
                                                                              self.z,
                                                                              self.nx,
                                                                              self.ny,
                                                                              self.nz)

parser = argparse.ArgumentParser()
parser.add_argument("starting_atom_index", help="the first atom index for the wall particles (e.g. one past the last atom index in the target file)")
parser.add_argument("atom_type", help="the atom type of the wall particles")
parser.add_argument("mol_id", help="the molecule id of the wall particles")
parser.add_argument("mass", help="the mass of each wall particle")
args = parser.parse_args()

atom_idx = int(args.starting_atom_index)
atom_type = int(args.atom_type)
mol_id = int(args.mol_id)
mass = float(args.mass)

atoms_expected = 0
line_idx = -1
#atoms = {} # indexed by atom id

for line in sys.stdin:
    line_idx += 1
    if line_idx == 0:
        # the first line contains the number of atoms in the xyz file
        atoms_expected = int(line)
    elif line_idx == 1:
        # the second line is reserved for comments: skip it
        continue
    else:
        atom = AtomBond(line, atom_idx, mol_id, atom_type)
        atom_idx += 1
        print(atom)
