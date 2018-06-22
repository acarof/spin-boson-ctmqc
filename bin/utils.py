#!/usr/bin/python
# -*- coding: latin-1 -*-

import numpy as np


def find_ctmqc_path():
    def create_ctmqc_path():
        print "The file ctmqc.path doesn't exist, please provide the path for ctmqc:"
        path = raw_input('> ')
        ctmqc_path = open('bin/ctmqc.path', 'w')
        ctmqc_path.write(path)
        ctmqc_path.close()
        return open('bin/ctmqc.path', 'r')
    try:
        ctmqc_path = open('bin/ctmqc.path', 'r')
    except:
        ctmqc_path = create_ctmqc_path()
    return ctmqc_path.readline()


def write_input(xpoints, friction, temperature, path, model, algorithm, final_time, dt, dump, ntraj):
    results = """
&SYSTEM
  MODEL_SYSTEM = '%s'
  N_DOF = 1           ! Tully s models are 1D
  X_POINTS = %s     ! check with wc -l *_bopes.dat
  Y_POINTS = 1
  Z_POINTS = 1
  NSTATES = 2
  EL_BASIS = 'adiabatic'
  DIA_TO_AD = 'y'
/
&DYNAMICS
  ALGORITHM = '%s' ! Ehrenfest, CTeMQC
  FINAL_TIME = %s
  DT = %s
  DUMP = %s
  INITIAL_BOSTATE = 0
  INITIAL_DIASTATE = 1
  NTRAJ = %s
  R_INIT = -10.0
  K_INIT = 25.0
  SIGMA_INIT = 0.0      ! if sigma = 0.0, then sigma = 20/k0 in the code
  MASS_INPUT = 1.0
  VISCOSITY = %s
  TEMPERATURE = %s
/
&EXTERNAL_FILES
  PATH_TO_POTENTIALS = "./spin_boson_surfaces_nacv/"
  POSITIONS_FILE = "./config_init/positions.dat"
  MOMENTA_FILE = "./config_init/velocities.dat"
/

    """ % ( model, xpoints, algorithm, final_time, dt, dump, ntraj,  friction, temperature,)
    with open("%s/input.in" % path, "w") as file_:
        file_.write(results)

def write_files(positions, mass, omega, epsilon_0, shift, coupling, temperature, path='.'):
    surf_ground = open("%s/1_bopes.dat" % path, 'w')
    surf_excited = open("%s/2_bopes.dat" % path, 'w')
    nacv_file = open("%s/nac1-12_x" % path, 'w')
    transfo_file = open("%s/transformation_matrix" % path, 'w')
    check_transo = open("%s/check_transfo.dat" % path, 'w')
    for position in positions:
        number = adiab_surfaces(position, mass, omega, epsilon_0, shift, coupling, True)
        surf_ground.write("%s  %s\n" % (number, position))
        number = adiab_surfaces(position, mass, omega, epsilon_0, shift, coupling, False)
        surf_excited.write("%s  %s\n" % (number, position))
        number = nacv(position, mass, omega, epsilon_0, shift, coupling)
        nacv_file.write("%s  %s\n" % (number, position))
        f1 = transfo_matrix(position, mass, omega, epsilon_0, shift, coupling, "f1")
        f2 = transfo_matrix(position, mass, omega, epsilon_0, shift, coupling, "f2")
        g1 = transfo_matrix(position, mass, omega, epsilon_0, shift, coupling, "g1")
        g2 = transfo_matrix(position, mass, omega, epsilon_0, shift, coupling, "g2")
        transfo_file.write(" %s  %s  %s  %s %s\n " % (f1, f2, g1, g2, position))
        check_transo.write("%s %s %s\n" % ( (f1**2+g1**2), (f2**2 + g2**2), f1*g1+f2*g2 ))
    surf_ground.close()
    surf_excited.close()
    nacv_file.close()
    transfo_file.close()
    check_transo.close()

def write_initial(mass, omega, epsilon_0, shift, coupling, temperature, path='.', seed=0):
    np.random.seed(seed)
    k_b = 0.00000316679085
    minimum = shift / (mass * omega ** 2)
    mean = -minimum
    sigma = np.sqrt(k_b*temperature / (mass * omega ** 2))
    distrib = list(np.random.normal(mean, sigma, 1000))
    with open("%s/config_init/positions.dat" % path, 'w') as file_:
        file_.write("%s" % '\n'.join(map(str,distrib)))

    mean = 0.0
    sigma = np.sqrt(k_b*temperature / (mass ))
    distrib = list(np.random.normal(mean, sigma, 1000))
    with open("%s/config_init/velocities.dat" % path, 'w') as file_:
        file_.write("%s" % '\n'.join(map(str,distrib)))

def adiab_surfaces(position, mass, omega, epsilon_0, shift, coupling, ground):
    if ground:
        sign = -1
    else:
        sign = 1
    ener = 0.5*mass*omega**2*position**2
    ener +=  - 0.5*epsilon_0
    root = np.sqrt((0.5*epsilon_0 + shift*position)**2 + coupling**2)
    ener += sign*root
    return ener

def nacv(position, mass, omega, epsilon_0, shift, coupling):
    num = 0.5*coupling*shift
    denom = (0.5*epsilon_0 + shift*position)**2 + coupling**2
    return num/denom


def transfo_matrix(position, mass, omega, epsilon_0, shift, coupling, indices):
    if indices == "f1":
        sign = (1,-1)
    elif indices == "f2":
        sign = (1,1)
    elif indices == "g1":
        sign = (-1,1)
    elif indices == "g2":
        sign = (1,-1)
    else:
        print "%s is not possible" % indices
        raise SystemError
    result = 0.5
    num = 0.5*(0.5*epsilon_0 + shift*position)
    denom = np.sqrt( (0.5*epsilon_0 + shift*position)**2 + coupling**2  ) 
    result += sign[1]*num/denom 
    result = sign[0]*np.sqrt(result)
    return result
