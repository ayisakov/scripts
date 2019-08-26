#!/usr/bin/env python3

# This script can be used to determine the density input to the Fortran-based
# chain creation tool (chain.f) that is distributed with LAMMPS in order to
# perform "sandwich" confinement simulations of polymer chains. The density
# input to chain.f must be higher than the desired density because the walls
# that comprise the confinement will need to be created at some distance away
# from the chains to aid in initial equilibration, thus reducing the density
# of the chains created by chain.f. Since the chain creation tool has no
# awareness of the walls, which will be created later in a separate step,
# this script becomes necessary. The extra space, or "padding" between the
# chains and the walls is necessary to make room for the particles that make
# up the walls. The wall particles are assumed to have a simple cubic unit
# cell with spacing = 0.8sigma

import argparse
import math
parser = argparse.ArgumentParser()
parser.add_argument("num_atoms", help="the number of chain atoms")
parser.add_argument("reduced_density", help="the desired reduced density of the chains only (L-J units, with sigma=1 assumed)")
parser.add_argument("padding", help="the extra amount of space on both sides of the chains in the direction normal to the confinement plane (h/2 on each side), units of sigma, with sigma=1 assumed")
args = parser.parse_args()

n = float(args.num_atoms)
rhostar = float(args.reduced_density)
padding = float(args.padding)

def alpha(n, rhostar, padding):
    one = 9. * padding**3. * n**2. * rhostar**4.
    radical = 27. * padding**6. * n**4. * rhostar**8. - 4. * padding**9. * n**3. * rhostar**9.
    two = math.sqrt(3.) * math.sqrt(radical)
    return one + two

def adjusted_density(n, rhostar, padding):
    al = alpha(n, rhostar, padding)**(1./3.)
    one = al * 2.**(-1./3.) * 3.**(-2./3.) * 1./n
    two = (2./3.)**(1./3.) * padding**3. * rhostar**3. * 1./al
    return one + two + rhostar

print(adjusted_density(n, rhostar, padding))
