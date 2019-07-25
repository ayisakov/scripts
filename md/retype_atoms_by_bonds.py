#!/usr/bin/env python3

# Usage:
# retype_atoms_by_bonds.py <type> <# of bonds> <new type> [<type> <# of bonds> <new type>][...] < data.in > data.out
# set the atom type to new type for atoms of original type with the specified number of bonds in a LAMMPS
# data file

import sys
from enum import Enum

# a Kremer-Grest "atom"
class AtomKG:
#    def __init__(self, _id, _mol_id, _type, _x, _y, _z, _nx, _ny, _nz):
#        self.atom_id = _id
#        self.mol_id = _mol_id
#        self.atom_type = _type
#        self.x = _x
#        self.y = _y
#        self.z = _z
#        self.nx = _nx
#        self.ny = _ny
#        self.nz = _nz
#        self.nbonds = 0
    def __init__(self, atom_line):
        contents = atom_line.split(' ')
        num_values = len(contents)
        if num_values != 9:
            raise Exception("An atom line is expected to have 9 values: {0:d} found.".format(num_values))
        self.atom_id = int(contents[0])
        self.mol_id = int(contents[1])
        self.atom_type = int(contents[2])
        self.x = float(contents[3])
        self.y = float(contents[4])
        self.z = float(contents[5])
        self.nx = int(contents[6])
        self.ny = int(contents[7])
        self.nz = int(contents[8])
        self.nbonds = 0

    def __str__(self):
        return "{0:d} {1:d} {2:d} {3:f} {4:f} {5:f} {6:d} {7:d} {8:d}".format(self.atom_id,
                                                                              self.mol_id,
                                                                              self.atom_type,
                                                                              self.x,
                                                                              self.y,
                                                                              self.z,
                                                                              self.nx,
                                                                              self.ny,
                                                                              self.nz)

class Bond:
    def __init__(self, bond_line):
        contents = bond_line.split(' ')
        num_values = len(contents)
        if num_values != 4:
            raise Exception("A bond line is expected to have 4 values: {0:d} found.".format(num_values))
        self.bond_id = int(contents[0])
        self.bond_type = int(contents[1])
        self.atom_1 = int(contents[2])
        self.atom_2 = int(contents[3])

    def __str__(self):
        return "{0:d} {1:d} {2:d} {3:d}".format(self.bond_id,
                                                self.bond_type,
                                                self.atom_1,
                                                self.atom_2)

class Section(Enum):
    preamble = 1
    atoms = 2
    other = 3
    bonds = 4
    rest = 5

found_atoms_section = False
found_atom_records = False

atoms = {} # indexed by atom id
bonds = [] # in order in which they were read
# output lines, indexed by section

data_out = { Section.preamble: [], Section.atoms: [], Section.other: [],
             Section.bonds: [], Section.rest: []}

def copy_line(line, _section):
    data_out[_section].append(line)

section = Section.preamble

for line in sys.stdin:
    if line == "": # EOF
        break

    if section == Section.preamble:
        # look for the Atoms section
        if line.split(' ')[0] != "Atoms":
            copy_line(line)
        else:
            section = Section.atoms

    else if section == Section.atoms:
        # populate atom records and look for the end of the section
        if line == "\n":
            if len(atoms) > 0:
                # end of atoms section
                section = Section.other
        else:
            # populate the atoms dictionary
            atom = AtomKG(line)
            atoms[atom.atom_id] = atom

    else if section == Section.other:
        # look for the Bonds section

    else if section == Section.bonds:
        # populate the bond records and look for the end of the section
        # retype the relevant atoms
        # print the atoms and bonds sections

    else:
        # copy over each line
        copy_line(line)

    if found_atoms_section:
        if found_atom_records:
            if line == "\n": # done with atom records
                found_atom_section = False

                continue
        else:
            # look for the first atom record
            if line == "\n": continue
            found_atom_records = True
    else:
        if line.split(' ')[0] != "Atoms":

