from __future__ import division
import numpy as np
import hoomd
from hoomd import md
from numpy import linalg as la
from Loggable import Loggable


class NetworkLJ(Loggable):

    def __init__(self, log_list=None, energy=15):
        super(NetworkLJ, self).__init__(log_list)

        self.log_values = ['pair_lj_energy']

        self.epsilon = [1, 1, 1, 1]
        self.sigma = [0.47, 0.47, 0.47, .47]
        self.names = ["BB", "A", "B", "FL"]

        self.lj_pair = None

        self.energy = energy

    def set_lj(self, neighbor_list, system):
        cut = 3 *.47
        self.lj_pair = hoomd.md.pair.lj(r_cut=cut, nlist=neighbor_list)

        self.add_to_logger()
        for t1 in system.particles.types:
            for t2 in system.particles.types:
                t1 = str(t1)
                t2 = str(t2)
                if t1 in self.names and t2 in self.names:
                    eps = np.sqrt(self.epsilon[self.names.index(t1)] * self.epsilon[self.names.index(t2)])
                    sig = (self.sigma[self.names.index(t1)] + self.sigma[self.names.index(t2)]) / 2
                    WCA_cut = 2 ** (1.0 / 6.0) * sig
                    if t1 == "A" and t2 == "A":
                        self.lj_pair.pair_coeff.set(str(t1), str(t2), epsilon=self.energy, sigma=.38 / 2 ** (1.0 / 6.0),
                                                    r_cut=cut)
                    else:
                        self.lj_pair.pair_coeff.set(str(t1), str(t2), epsilon=eps, sigma=sig, r_cut=WCA_cut)
                else:
                    self.lj_pair.pair_coeff.set(str(t1), str(t2), epsilon=0, sigma=1, r_cut=1)

    def disable(self):
        self.lj_pair.disable()

    def is_center(self, string):

        if len(string) >= 6 and string[:6] == 'center':
            return True
        return False


class LJ(Loggable):
    def __init__(self, log_list=None, energy=1):

        super(LJ, self).__init__(log_list)

        self.log_values = ['pair_lj_energy']

        self.epsilon = [0.5, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        #self.epsilon =  [1,   1, 1, 1, 1, 1 ,1]
        #self.sigma = [1.0, 0.5, 0.5, 0.5, 1.0, 0.4, 0.4]
        self.sigma = [1.0, 1.0, 1.0, 1.0, 1.0, 0.6, 0.6, 0.6, 0.6, 0.6]
        #self.sigma = [0.5, 0.5, 0.5, 0.5, 0.5, 0.3, 0.3]
        self.names = ['Bb', 'B', 'L', 'QM', 'QP', 'QMsmall', 'QMi', 'QPi', 'Q2Pi', 'Q2Mi']

        self.lj_pair = None

        self.energy = energy

    def set_lj(self, neighbor_list, system):
        cut = 3
        self.lj_pair = hoomd.md.pair.lj(r_cut=cut, nlist=neighbor_list)
        self.add_to_logger()
        for t1 in system.particles.types:
            for t2 in system.particles.types:
                t1 = str(t1)
                t2 = str(t2)
                if t1[0] == 'B' and t2[0] == 'B':
                    if len(t1) > 1 and len(t2)>1:
                        self.lj_pair.pair_coeff.set(str(t1), str(t2), epsilon=1, sigma=1, r_cut=3)
                    else:
                        sig = (self.sigma[self.names.index(t1)] + self.sigma[self.names.index(t2)]) / 2
                        cut = 3 * sig
                        eps = np.sqrt(self.epsilon[self.names.index(t1)] * self.epsilon[self.names.index(t2)])
                        self.lj_pair.pair_coeff.set(str(t1), str(t2), epsilon=self.energy * eps, sigma=sig, r_cut=cut)
                elif t1 in self.names and t2 in self.names:
                    eps = np.sqrt(self.epsilon[self.names.index(t1)] * self.epsilon[self.names.index(t2)])
                    sig = (self.sigma[self.names.index(t1)] + self.sigma[self.names.index(t2)]) / 2
                    WCA_cut = 2 ** (1.0 / 6.0) * sig
                    self.lj_pair.pair_coeff.set(str(t1), str(t2), epsilon=eps, sigma=sig, r_cut=WCA_cut)
                else:
                    self.lj_pair.pair_coeff.set(str(t1), str(t2), epsilon=0, sigma=sig, r_cut=1)

    def disable(self):
        self.lj_pair.disable()

    def is_center(self, string):

        if len(string) >= 6 and string[:6] == 'center':
            return True
        return False


class ThreePM(Loggable):
    def __init__(self, log_list=None, eps_relative=None, lb=None):

        super(ThreePM, self).__init__(log_list)

        if eps_relative is None and lb is None:
            eps_relative = 80

        if eps_relative is None and lb is not None:
            eps_relative = 56/lb

        if lb is None and eps_relative is not None:
            lb = 56/ eps_relative
            #lb= 112/eps_relative


        #if lb is not None and eps_relative is not None and eps_relative * lb != 56:
         #   raise ValueError("lb and eps_relative are inconsistent. Must multiply to 56")


        self.log_values = ['pppm_energy']
        self.charge = [0, 0, 0, -1, 1, -1, -1, 1, 2, -2, 0]
        self.names = ['Bb', 'B', 'L', 'QM', 'QP', 'QMsmall', 'QMi', 'QPi', 'Q2Pi', 'Q2Mi',  'center']
        self.charge_adjusted = False
        # self.charge_unit = 1.6 * 10^-19 / sqrt( 4 pi eps_0 eps_r  4.11 * 10^-21 J/kt 10^-9 m )
        self.charge_unit = .837 * np.sqrt(80) / np.sqrt(eps_relative)

        self.object = None

    def set_charges(self, neighbor_list, system, lamda=1):

        found_charge = False
        for part in system.particles:
            if self.charge[self.names.index(part.type)] != 0:
                found_charge = True
                if not self.charge_adjusted:
                    part.charge = part.charge * self.charge_unit * lamda

        if found_charge:
            self.object = hoomd.md.charge.pppm(group=hoomd.group.charged(), nlist=neighbor_list)
            #self.object = hoomd.md.charge.pppm(group=hoomd.group.all(), nlist=neighbor_list)
            self.add_to_logger()
            #self.object.set_params(Nx=32, Ny=32, Nz=32, order=5, rcut=5 * 2 ** (1/6))
            self.object.set_params(Nx=64, Ny=64, Nz=64, order=6, rcut=3 * 2 ** (1 / 6))
            self.charge_adjusted = True

        if self.charge_adjusted:
            self.charge_unit = 1

    def turn_off_charge(self):

        if self.object is not None:
            self.object.disable()
            self.remove_from_logger()

    def turn_on_charge(self, lamda=1):

        if self.object is not None:
            self.object.enable()
            self.add_to_logger()


def LJRepulsive_pair(r,rmin, rmax, sigma, epsilon):

    if r < sigma:
        V = epsilon * ((sigma / r) ** 12 - 1)
        F = epsilon * 12 * (sigma/r) ** 13
    else:
        V = 0
        F = 0
    return (V,F)


def quad(r, sigma, epsilon):
    if r < (sigma + .25) and r > (sigma - .25):
        V = 16 * epsilon * (r - 1) ** 2 - epsilon
        F = - (16 * epsilon * (2 * r - 2))
    else:
        V = 0
        F = 0
    return (V, F)


def H_potential(r, rmin, rmax, sigma, epsilon):
    if r > .5:
        return quad(r, sigma, epsilon)
    else:
        return LJRepulsive_pair(r,0,0, .5, 1)


def gaussian(std, mean, x):

    return np.exp(-(x-mean)**2 / (2 * std**2))

def gaussian_prime(std, mean, x):

    return 2 * (x-mean)/(2 * std**2) * gaussian(std, mean, x)

def double_well(r, rmin, rmax, well1, well2):

    if r > 1:
        V =0
        F =0
    elif r > .65:
        V = -well1 * gaussian(.05, .8,  r)
        F = - well1 * gaussian_prime(.05, .8, r)

    elif r > .35:
        V = -well2 * gaussian(.05, .5, r)
        F = - well2 * gaussian_prime(.05, .5, r)
    else:
        return LJRepulsive_pair(r, 0, 0, .35, 1)
    return (V,F)

def double_well2(r, rmin, rmax, well1, well2, barrier):

    if r > 1:
        V =0
        F =0
    elif r > .72:
        V = -well1 * gaussian(.04, .84,  r)
        F = - well1 * gaussian_prime(.04, .84, r)
    elif r > .48:
        V = barrier * gaussian(.04, .6, r)
        F = barrier * gaussian_prime(.04, .6, r)
    elif r > .24:
        V = -well2 * gaussian(.04, .36, r)
        F = - well2 * gaussian_prime(.04, .36, r)
    else:
        return LJRepulsive_pair(r, 0, 0, .24, 1)
    return (V,F)

def ci_release(r, rmin, rmax, sphere_radius, effective_charge, debye_length):

    V = 0
    F = 0
    if r > sphere_radius + 1:
        V = - 2 * (1 - effective_charge)  * np.exp(-(r - (sphere_radius + 1))/ debye_length)
        F = - 2 * (1 - effective_charge) / debye_length * np.exp(-(r - (sphere_radius + 1))/ debye_length)
    else:
        V = - 2 * (1 - effective_charge)
        F = 0

    print(r,V,F, sphere_radius)
    return V, F

def lj_sphere(r, rmin, rmax, sphere_radius, effective_charge, debye_length):

    V = 0
    F = 0
    if r > sphere_radius + 1:
        V = - 1 + np.power((r - 1 - sphere_radius), 6)
        F = -6 * np.power((r - 1 - sphere_radius), 5)
    else:
        V = - 1
        F = 0

    return V, F