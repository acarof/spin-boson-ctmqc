from utils import *
import numpy as np
from itertools import product as iterprod
import os

mass = 1
omega = 3.5E-4
epsilon_0 = 1.5E-2
reorga = 2.39E-2
shift = np.sqrt( 0.5*reorga*mass*omega**2 )
coupling = 1.49E-5
temperature = 9.5E-4

list_param = [
    ['MASS', 1],
    ['OMEGA', 3.5E-4, 1000],
    ['EPSILON_0', 1.5E-2],
    ['REORGA', 2.39E-2],
   # ['COUPLING', 1.49E-5],
    ['COUPLING', 1.49E-5, 1, 2, 3],
    ['TEMPERATURE', 9.5E-4],
]

mega_list = []
second_list = [sublist[1:] for sublist in list_param]
total_list = list(iterprod(*second_list))
for sublist in total_list:
    subdict = {}
    for index in range(len(sublist)):
        subdict.update({
            list_param[index][0]: sublist[index]
        })
    mega_list.append(subdict)
print len(mega_list)


for index, dict_ in enumerate(mega_list):
    mass = dict_['MASS']
    omega = dict_['OMEGA']
    epsilon_0 = dict_['EPSILON_0']
    reorga = dict_['REORGA']
    coupling = dict_['COUPLING']
    temperature = dict_['TEMPERATURE']

    shift = np.sqrt(0.5 * reorga * mass * omega ** 2)
    minimum = shift/(mass*omega**2)
    positions = np.arange(-2*minimum,2*minimum,0.1)

    name_dir = 'run-ctmqc-%s' % index
    os.mkdir(name_dir)
    write_files(positions, mass, omega, epsilon_0, shift, coupling, temperature, path=name_dir)
    write_initial(mass, omega, epsilon_0, shift, coupling, temperature, path=name_dir)