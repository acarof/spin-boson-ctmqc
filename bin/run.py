#!/usr/bin/python
# -*- coding: latin-1 -*-

from utils import *
import numpy as np
from itertools import product as iterprod
import os, time, sys, random
import subprocess, hashlib
from list_param import *

ctmqc = find_ctmqc_path()




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


def get_md5name():
    complete_time = time.strftime("%y%m%d%H%M%S", time.localtime())
    md5 = hashlib.md5()
    md5.update(complete_time + "%s" % random.randint(1, 1E9))
    md5name = md5.hexdigest()
    return md5name


def run_ctmqc(dict_):
    mass = dict_['MASS']
    omega = dict_['OMEGA']
    epsilon_0 = dict_['EPSILON_0']
    reorga = dict_['REORGA']
    coupling = dict_['COUPLING']
    temperature = dict_['TEMPERATURE']
    friction = dict_['VISCOSITY']
    model = dict_['MODEL']
    algorithm = dict_['ALGORITHM']
    final_time = dict_['FINAL_TIME']
    dt = dict_['DT']
    dump = dict_['DUMP']
    ntraj  = dict_['NTRAJ']
    grid = dict_.get('GRID', 0.1)

    shift = np.sqrt(0.5 * reorga * mass * omega ** 2)
    minimum = shift/(mass*omega**2)
    positions = np.arange(-4*minimum,4*minimum, grid)

    name_dir = 'run-ctmqc-%s' % get_md5name()
    os.mkdir(name_dir)

    for dir in ['output', 'output/histo', 'output/coeff', 'output/trajectories', 'output/density',
                'spin_boson_surfaces_nacv', 'config_init']:
        os.mkdir('%s/%s' % (name_dir, dir))

    write_files(positions, mass, omega, epsilon_0, shift, coupling, temperature,
                path='%s/spin_boson_surfaces_nacv/' % name_dir)
    write_initial(mass, omega, epsilon_0, shift, coupling, temperature, path=name_dir)
    write_input(xpoints=len(positions), friction=friction, temperature=temperature, path=name_dir,
                model=model, algorithm=algorithm, final_time=final_time, dt=dt, dump=dump, ntraj=ntraj)

    with open("%s/parameters.dat" % name_dir, 'w') as file_:
        for key, value in dict_.iteritems():
            file_.write("%s  %s\n" % (key, value))

    os.chdir(name_dir)

    complete_time = time.strftime("%y%m%d%H%M%S", time.localtime())
    # This section sets up a list with the aprun command
    runcommand = []
    print "CTMQC STARTS AT: " + complete_time
    execName = ctmqc
    #inputName = "run.inp"
    #outputName = "run.log"
    # Add executable name to the command
    runcommand.append(execName)
    #runcommand += ['-i', inputName]
    print runcommand
    #logFile = open(outputName, 'w')
    inFile = open("input.in")
    stderr = subprocess.call(runcommand, stdin=inFile,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    os.chdir('..')

# RUN THE CALCULATIONS, SERIE OR PARALLEL ACCORDING TO THE NWORKER VARIABLE
try:
    nworker = int(sys.argv[1])
except:
    nworker = 1

print "nworker is %s" % nworker

if nworker == 1:
    print 'Use serial'
    for dict_ in mega_list:
        run_ctmqc(dict_)
elif nworker == 0:
    from multiprocessing import Pool, cpu_count
    pool = Pool(cpu_count())
    print 'Use parallel with %s processors' % cpu_count()
    pool.map(run_ctmqc, mega_list)
elif nworker > 1:
    from multiprocessing import Pool
    pool = Pool(nworker)
    print 'Use parallel with %s processors' % nworker
    pool.map(run_ctmqc, mega_list)





