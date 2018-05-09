import numpy as np


def write_files(positions, mass, omega, epsilon_0, shift, coupling, temperature, path='.'):
    surf_ground = open("%s/1_bopes.dat" % path, 'w')
    surf_excited = open("%s/2_bopes.dat" % path, 'w')
    nacv_file = open("%s/nac1_12_x" % path, 'w')
    transfo_file = open("%s/transformation_matrix.dat" % path, 'w')
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

def write_initial(mass, omega, epsilon_0, shift, coupling, temperature, path='.'):
    minimum = shift / (mass * omega ** 2)
    mean = -minimum
    sigma = np.sqrt(temperature / (mass * omega ** 2))
    distrib = list(np.random.normal(mean, sigma, 1000))
    with open("%s/initial_positions.dat" % path, 'w') as file_:
        file_.write("%s" % '\n'.join(map(str,distrib)))

    mean = 0.0
    sigma = np.sqrt(temperature / (mass ))
    distrib = list(np.random.normal(mean, sigma, 1000))
    with open("%s/initial_velocities.dat" % path, 'w') as file_:
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
