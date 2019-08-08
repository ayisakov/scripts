#!/usr/bin/env python3

# Usage:
# retype_atoms_by_bonds.py <type> <# of bonds> <new type> [<type> <# of bonds> <new type>][...] < data.in > data.out
# set the atom type to new type for atoms of original type with the specified number of bonds in a LAMMPS
# data file

import sys
from enum import Enum

# TODO: get these parameters from the arguments list
# atom_type: (nbonds, new_type)
retype = { 1: (1, 3), 2: (1, 4)}

# TODO: add masses for the new atom types

# a Kremer-Grest "atom"
class AtomKG:
    def __init__(self, atom_line):
        contents = atom_line.split()
        num_values = len(contents)
        if num_values != 9:
            raise Exception("An atom line is expected to have 9 values: {0:d} found." \
            "\nLine: {1:s}".format(num_values, atom_line))
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
        contents = bond_line.split()
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


data_out = { Section.preamble: [], Section.atoms: ["Atoms # bond", ""], Section.other: [],
             Section.bonds: ["Bonds", ""], Section.rest: []}

def copy_line(line, _section):
    data_out[_section].append(line)

def find_section(line, section):
    if line[:(len(section))] == section:
#    if line.split()[0] == section:
        return True
    return False

section = Section.preamble

retype_done = False
def do_retype():
    # retype the relevant atoms
    for k in atoms.keys():
        atom = atoms[k]
        if atom.atom_type in retype:
            #print("DEBUG: atom found in retype")
            r = retype[atom.atom_type]
            if atom.nbonds == r[0]:
                #print("DEBUG: changing atom type")
                atom.atom_type = r[1]
                update_n_atom_types(atom.atom_type)
    # store the atoms and bonds sections
    for atom in atoms.values():
        data_out[Section.atoms].append(str(atom))
    data_out[Section.atoms].append("")
    #print("DEBUG: {0:d} bonds".format(len(bonds)))
    for bond in bonds:
        data_out[Section.bonds].append(str(bond))
    retype_done = True

n_types = 0
def update_n_atom_types(n):
    global n_types
    if n > n_types:
        n_types = n

for line in sys.stdin:
    # TODO: this check looks unnecessary
    if line == "": # EOF
        break

    if section == Section.preamble:
        # look for the Atoms section
        if not find_section(line, "Atoms"):
            copy_line(line, section)
        else:
#            #print("DEBUG: Found section atoms")
            section = Section.atoms
            continue

    elif section == Section.atoms:
        # populate atom records and look for the end of the section
        if line == "\n":
#            #print("DEBUG: Found empty line in section 'Atoms'")
            if len(atoms) > 0:
                # end of atoms section
                section = Section.other
        else:
            # populate the atoms dictionary
#            #print("DEBUG: Found an atom")
            atom = AtomKG(line)
            atoms[atom.atom_id] = atom
            update_n_atom_types(atom.atom_type)

    elif section == Section.other:
#        #print("DEBUG: in section 'other'")
        # look for the Bonds section
        if not find_section(line, "Bonds"):
            copy_line(line, section)
            #print("DEBUG: line is: "+line)
        else:
            #print("DEBUG: Found section 'Bonds'")
            section = Section.bonds
            continue

    elif section == Section.bonds:
        # populate the bond records and look for the end of the section
        #print("DEBUG: in section 'Bonds'; line is: "+line)
        if line == "\n":
            if len(bonds) > 0:
                #print("DEBUG: end of section 'Bonds'")
                # end of bonds section
                section = Section.rest
                do_retype()
        else:
            #print("DEBUG: Found a bond")
            bond = Bond(line)
            bonds.append(bond)
            atoms[bond.atom_1].nbonds += 1
            atoms[bond.atom_2].nbonds += 1
    else:
        # copy over each line
        copy_line(line, section)

if not retype_done:
    do_retype()

# Update the number of atom types in the preamble
preamble = data_out[Section.preamble]
for i in range(0, len(preamble) - 1):
    l = preamble[i]
    if "atom types" in l:
        preamble[i] = "{0:d} atom types".format(n_types)

# Print the sections in order
for s in [Section.preamble, Section.atoms, Section.other, Section.bonds, Section.rest]:
    for line in data_out[s]:
        print(line.rstrip())
