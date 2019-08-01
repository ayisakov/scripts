#!/usr/bin/env python3

# Usage:
# unwrap_periodic.py < data.in > data.out
# Unwraps all atom coordinates from a LAMMPS data file such that all coordinates
# with image flags are replaced with equivalent coordinates in image 0. The box
# is resized to accommodate the new coordinates.

import sys
import copy
from enum import Enum
import numpy as np

# a Kremer-Grest "atom"
class AtomKG:
    def __init__(self, atom_line):
        self.pos = [0., 0., 0.]
        self.image = [0, 0, 0]
        contents = atom_line.split(' ')
        num_values = len(contents)
        if num_values != 9:
            raise Exception("An atom line is expected to have 9 values: {0:d} found.\nLine was: {1:s}".format(num_values, atom_line))
        self.atom_id = int(contents[0])
        self.mol_id = int(contents[1])
        self.atom_type = int(contents[2])
        self.pos[0] = float(contents[3])
        self.pos[1] = float(contents[4])
        self.pos[2] = float(contents[5])
        self.image[0] = int(contents[6])
        self.image[1] = int(contents[7])
        self.image[2] = int(contents[8])
        self.nbonds = 0

    def __str__(self):
        return "{0:d} {1:d} {2:d} {3:f} {4:f} {5:f} {6:d} {7:d} {8:d}".format(self.atom_id,
                                                                              self.mol_id,
                                                                              self.atom_type,
                                                                              self.pos[0],
                                                                              self.pos[1],
                                                                              self.pos[2],
                                                                              self.image[0],
                                                                              self.image[1],
                                                                              self.image[2])

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
        self.length = None

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

class Box:
    def __init__(self, xlo, xhi, ylo, yhi, zlo, zhi):
        self.dim = [[0., 0.], [0., 0.], [0., 0.]]
        self.dim[0][0] = xlo
        self.dim[0][1] = xhi
        self.dim[1][0] = ylo
        self.dim[1][1] = yhi
        self.dim[2][0] = zlo
        self.dim[2][1] = zhi
    def xlen(self):
        return self.length(0)
    def ylen(self):
        return self.length(1)
    def zlen(self):
        return self.length(2)
    def length(self, dimension):
        return self.dim[dimension][1] - self.dim[dimension][0]
    def __str__(self):
        return "{0:9f} {1:9f} xlo xhi\n{2:9f} {3:9f} ylo yhi\n{4:9f} {5:9f} zlo zhi".format(self.dim[0][0],
                                                                                            self.dim[0][1],
                                                                                            self.dim[1][0],
                                                                                            self.dim[1][1],
                                                                                            self.dim[2][0],
                                                                                            self.dim[2][1])

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
#    if line.split(' ')[0] == section:
        return True
    return False

section = Section.preamble

def fullpos(box, atom):
    def fulldim(axis, reldim, imagenum):
#        print("axis = ", axis, "reldim = ", reldim, "imagenum = ", imagenum)
        return reldim + box.length(axis) * imagenum
    return [fulldim(ax, atom.pos[ax], atom.image[ax]) for ax in [0, 1, 2]]

def calculate_distance(box, atom1, atom2):
    pos1 = np.array(fullpos(box, atom1))
#    print('Debug: atom {0:d} (1) y-position = {1:.3f}'.format(atom1.atom_id, pos1[1]))
    pos2 = np.array(fullpos(box, atom2))
#    print('Debug: atom {0:d} (2) y-position = {1:.3f}'.format(atom2.atom_id, pos2[1]))
    dist = pos2 - pos1
#    print('Debug: y-distance = {0:.3f}'.format(dist[1]))
    return np.sqrt(np.sum(np.square(dist)))


n_types = 0

thebox = Box(0., 0., 0., 0., 0., 0.)
bond_length_min = 99999999.
bond_length_max = 0.

newbox = Box(0., 0., 0., 0., 0., 0.)

def unwrap_atom_coordinates():
    padding = 0.1
    for atom_id in atoms.keys():
        atom = atoms[atom_id]
        atom.pos = fullpos(thebox, atom)
        atom.image = [0, 0, 0]
        # resize the box to fit all coordinates
        for ax in [0, 1, 2]:
            if atom.pos[ax] <= newbox.dim[ax][0]:
                newbox.dim[ax][0] = atom.pos[ax] - padding
            if atom.pos[ax] >= newbox.dim[ax][1]:
                newbox.dim[ax][1] = atom.pos[ax] + padding
        data_out[Section.atoms].append(str(atom))
    data_out[Section.atoms].append("")

def process_bonds():
    for bond in bonds:
        data_out[Section.bonds].append(str(bond))

distances_done = False

for line in sys.stdin:
    # TODO: this check looks unnecessary
    if line == "": # EOF
        break

    if section == Section.preamble:
        # look for the Atoms section
        if not find_section(line, "Atoms"):
            # look for the box dimensions
            linesp = line.split(' ')
            if "xhi" in line:
#                print("Debug: found x box line")
                thebox.dim[0][0] = float(linesp[0])
                thebox.dim[0][1] = float(linesp[1])
            elif "yhi" in line:
#                print("Debug: found y box line")
                thebox.dim[1][0] = float(linesp[0])
                thebox.dim[1][1]= float(linesp[1])
            elif "zhi" in line:
#                print("Debug: found z box line")
                thebox.dim[2][0] = float(linesp[0])
                thebox.dim[2][1] = float(linesp[1])
            else:
                # Do not copy the box lines -- will write our own
                copy_line(line, section)

        else:
#            #print("DEBUG: Found section atoms")
            newbox = copy.deepcopy(thebox)
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
                unwrap_atom_coordinates()
                process_bonds()
        else:
            #print("DEBUG: Found a bond")
            bond = Bond(line)
            bonds.append(bond)
            atoms[bond.atom_1].nbonds += 1
            atoms[bond.atom_2].nbonds += 1
    else:
        # copy over each line
        copy_line(line, section)

if not distances_done:
    unwrap_atom_coordinates()
    process_bonds()

# Update the number of atom types in the preamble
#preamble = data_out[Section.preamble]
#for i in range(0, len(preamble) - 1):
#    l = preamble[i]
#    if "atom types" in l:
#        preamble[i] = "{0:d} atom types".format(n_types)

# Update preamble with new box information
idx_masses = 0
for line in data_out[Section.preamble]:
    if "Masses" in line:
        break
    idx_masses += 1
data_out[Section.preamble].insert(idx_masses, "")
data_out[Section.preamble].insert(idx_masses - 1, str(newbox))

# Print the sections in order
for s in [Section.preamble, Section.atoms, Section.other, Section.bonds, Section.rest]:
    for line in data_out[s]:
        print(line.rstrip())
# Print bond information
#for bond in bonds:
#    print(str(bond) + " length: " + str(bond.length))
#    if bond.length > 2.5:
#        print("WARNING: long bond; atom lines below: ")
#        print(atoms[bond.atom_1])
#        print(atoms[bond.atom_2])
#        print('')
#print("Shortest bond: {0:.3f}\nLongest bond: {1:.3f}".format(bond_length_min,
#                                                             bond_length_max))
#print("Box dimensions (x, y, z): {0:.2f}, {1:.2f}, {2:.2f}".format(thebox.xlen(),
#                                                                   thebox.ylen(),
#                                                                   thebox.zlen()))
