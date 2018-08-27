import sys
import cantera as ct
import numpy as np


cti = ct.Solution('boivin.cti')

fuel_species = 'H2'

m = cti.n_species

ifuel = cti.species_index(fuel_species)
io2 = cti.species_index('O2')

#print('THERMO PROPERTIES------------\n')
#phi = input('Enter stoichiometric ratio: ')
#phi = float(phi)

x = np.zeros(m,'d')

#x[ifuel] = phi
#stoich_o2 = cti.n_atoms(fuel_species,'H')
#x[io2] = stoich_o2

#Ti = input('Enter temperature: ')
#Ti = float(Ti)

cti.TPX = 1200,ct.one_atm, 'H2:2, O2:1'

#Creating the reactor
r = ct.IdealGasReactor(cti)
#Need to specify const pressure or const volume- we take constant pressure

#Create a reactor network

sim = ct.ReactorNet([r])
time = 0.0
states = ct.SolutionArray(cti, extra=['t'])

print('%10s %10s %10s %14s %14s' % ('t [s]','T [K]','P [Pa]','u [J/kg]','mass(kg)'))
for n in range(1000):
    time += 1.e-7
    sim.advance(time)
    states.append(r.thermo.state, t=time*1e3)
    print('%10.3e %10.3f %10.3f %14.6e %14.6e' % (sim.time, r.T,
                                           r.thermo.P, r.thermo.u,r.mass))



if '--plot' in sys.argv[1:]:
    import matplotlib.pyplot as plt
    plt.clf()
    plt.subplot(2, 2, 1)
    plt.plot(states.t, states.T)
    plt.xlabel('Time (ms)')
    plt.ylabel('Temperature (K)')
    plt.subplot(2, 2, 2)
    plt.plot(states.t, states.X[:,cti.species_index('OH')])
    plt.xlabel('Time (ms)')
    plt.ylabel('OH Mole Fraction')
    plt.subplot(2, 2, 3)
    plt.plot(states.t, states.X[:,cti.species_index('O2')])
    plt.xlabel('Time (ms)')
    plt.ylabel('O2 Mole Fraction')
    plt.subplot(2, 2, 4)
    plt.plot(states.t, states.X[:,cti.species_index('H2')])
    plt.xlabel('Time (ms)')
    plt.ylabel('H2 Mole Fraction')
    plt.tight_layout()
    plt.show()
else:
    print("To view a plot of these results, run this script with the option --plot")

