#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script makes DFT calculations to get free energy for the oxygen evolution reaction
DFT calculations are performed with GPAW and ASE is used for system design

It uses very simple models with several considerations to reduce the computational cost
Therefore, it is not useful for production and the energy values would be not very accurate

This is only created to describe how to use DFT calculations
for plotting free energy diagrams of an electrochemical reaction

You can read the first part of the tutorial here:
http://www.dyerga.org/blog/2019/02/09/practical-introduction-to-dft-for-electrocatalysis-1-free-energy-diagrams/

author: Daniel Martin-Yerga
email: dyerga@gmail.com

"""
from ase import Atoms
from ase.build import fcc111, molecule
from ase.constraints import FixAtoms
from ase.visualize import view
from ase.optimize import QuasiNewton
from gpaw import GPAW, PW, FermiDirac
import time


"""
Calculations to obtain the free energy of the Ni(111) slab
2.0 A are used for the vacuum, which is a low value, but higher values increased a lot the calculation time
"""

lattice_constant_Ni = 3.52
slab = fcc111("Ni", a=lattice_constant_Ni, size=(2,2,2))
slab.center(vacuum=2.0)
view(slab)
calculator = GPAW(xc='RPBE', mode=PW(350), kpts={'size': (4,4,1), 'gamma': True}, h=0.2, occupations=FermiDirac(0.1))
slab.set_calculator(calculator)
"""Here we fix the position of 7 of 8 atoms of the slab to facilitate the calculation"""
constraint = FixAtoms(indices=[0,1,2,3,5,6,7])
slab.set_constraint(constraint)


# We use the QuasiNewton algorithm for the structure optimization
# We save the optimized structure in the trajectory file
# Optimization will stop when force between atoms is lower than fmax: 0.05 eV
# We can also to record the time taken for the calculation
dyn = QuasiNewton(slab, trajectory='Ni.traj')
t = time.time()
dyn.run(fmax=0.05)
print('Calculation time: {} min.'.format((time.time() - t) / 60))

# Get tje final energy of the optimized slab structure
e_slab = slab.get_potential_energy()
print("Ni(111) energy: ", e_slab, " eV")



"""
Calculations to obtain the free energy of the OH adsorbed on the Ni slab
"""

# Creating the OH molecule and a new Ni slab
oh_molecule = Atoms('OH', positions=[(0, 0, 0), (0, -0.763, 0.596)])
slabNi = fcc111("Ni", a=lattice_constant_Ni, size=(2,2,2))

# Placing the molecule close to a top Ni atom, where it would bind
p = slabNi.positions[4]
oh_molecule.translate(p + (0, 0, 1.5))
slabOH = slabNi + oh_molecule
slabOH.center(vacuum=2.0)

# We fix again some atoms of the slab to speed up calculations
constraint = FixAtoms(indices=[0,1,2,3,5,6,7])
slabOH.set_constraint(constraint)
view(slabOH)

slabOH.set_calculator(calculator)

dynOH = QuasiNewton(slabOH, trajectory="OH_Ni.traj", )
t = time.time()
dynOH.run(fmax=0.05)
print('Calculation time: {} min.'.format((time.time() - t) / 60))

e_slabOH = slabOH.get_potential_energy()
print("OH on Ni(111) energy: ", e_slabOH, " eV")
print("bond Ni-O: ", slabOH.get_distance(4,8))
print("bond O-H: ", slabOH.get_distance(8,9))


"""
Calculations to get the free energy of the O intermediate adsorbed on the Ni surface
"""
o_molecule = Atoms('O', positions=[(0, 0, 0)])
slabNi = fcc111("Ni", a=lattice_constant_Ni, size=(2,2,2))
p = slabNi.positions[4]
o_molecule.translate(p + (0, 0, 1.5))
slabO = slabNi + o_molecule
slabO.center(vacuum=2.0)
constraint = FixAtoms(indices=[0,1,2,3,5,6,7])
slabO.set_constraint(constraint)
view(slabO)

slabO.set_calculator(calculator)
dynO = QuasiNewton(slabO, trajectory="O_Ni.traj", )
t = time.time()
dynO.run(fmax=0.05)
print('Calculation time: {} min.'.format((time.time() - t) / 60))

e_slabO = slabO.get_potential_energy()
print("O on Ni(111) energy: ", e_slabO, " eV")
view(slabO)
print("bond Ni-O: ", slabO.get_distance(4,8))


"""
Calculations to get the free energy of the OOH intermediate adsorbed on the Ni surface
"""
ooh_molecule = Atoms('OOH', positions=[(0, 0, 0), (0, 0, 1.4), (0, -0.763, 2.0)])
slabNi = fcc111("Ni", a=lattice_constant_Ni, size=(2,2,2))
p = slabNi.positions[4]
ooh_molecule.translate(p + (0, 0, 1.5))
slabOOH = slabNi + ooh_molecule
slabOOH.center(vacuum=2.0)
constraint = FixAtoms(indices=[0,1,2,3,5,6,7])
slabOOH.set_constraint(constraint)
view(slabOOH)

slabOOH.set_calculator(calculator)
dynOOH = QuasiNewton(slabOOH, trajectory="OOH_Ni.traj", )
t = time.time()
dynOOH.run(fmax=0.05)
print('Calculation time: {} min.'.format((time.time() - t) / 60))

e_slabO = slabOOH.get_potential_energy()
print("OOH on Ni(111) energy: ", e_slabO, " eV")
view(slabOOH)
print("bond Ni-O: ", slabOOH.get_distance(4,8))




"""
Calculations to get the free energy of H2O and H2 needed to calculate the free energies of the reaction
"""

calculator2 = GPAW(xc='RPBE', mode=PW(350), h=0.2, occupations=FermiDirac(0.1))
h2o_molecule = molecule("H2O")
h2o_molecule.set_cell(slab.get_cell())
h2o_molecule.center(vacuum=2.0)
h2o_molecule.set_calculator(calculator2)
dynH2O = QuasiNewton(h2o_molecule, trajectory='H2O.traj')
t = time.time()
dynH2O.run(fmax=0.05)
print('Calculation time: {} min.'.format((time.time() - t) / 60))
e_h2o = h2o_molecule.get_potential_energy()
print("H2O energy: ", e_h2o, " eV")
print("O-H bond lenght: ", h2o_molecule.get_distance(0,1))


h2_molecule = molecule("H2")
h2_molecule.set_cell(slab.get_cell())
h2_molecule.center(vacuum=2.0)
h2_molecule.set_calculator(calculator2)
dynH2 = QuasiNewton(h2_molecule, trajectory='H2.traj')
t = time.time()
dynH2.run(fmax=0.05)
print('Calculation time: {} min.'.format((time.time() - t) / 60))
e_h2 = h2_molecule.get_potential_energy()
print("H2 energy: ", e_h2, " eV")
print("H-H bond lenght: ", h2_molecule.get_distance(0,1))




