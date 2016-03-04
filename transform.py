'''
Code to do transformation from offset in 1 filter system to offset in another

Doing with mindset (currently) of:
given starspot modulation from Kepler, a spot delta T, and a stellar T_eff
 predict the ugrizJHK response
'''


import numpy as np
import matplotlib.pyplot as plt
from astropy.analytic_functions import blackbody_lambda
import pysynphot


def measure_ampl(flux, Tstart=5777.0, dTspot=500.0):
    '''
    Given a Kepler light curve (flux values in fractional flux?), assume stellar & spot temperatures,
    determine the filling factor of the spot relative to the stellar surface

    '''

    # input PHOENIX models

    # input Kepler filter
    Fdir = 'filters/'
    Ffile = 'kepler_response_hires1'

    # convolve spot & star models w/ filter

    # convert Delta Flux values in to linear combination of convolved models

    fillingfactor = np.zeros_like(flux)

    return fillingfactor


def transform(ffactor, Tstar=5777.0, dTspot=500.):
    '''
    Given (spot filling factor light curve, stellar & spot temperatures),
    compute the change in flux in all (ugriz) bands

    '''

    Fdir = 'filters/'
    Ffiles = {'u':Fdir + 'u.dat.txt', 'g':Fdir + 'g.dat.txt', 'r':Fdir + 'r.dat.txt',
              'i':Fdir + 'i.dat.txt', 'z':Fdir + 'z.dat.txt',
              'j':Fdir + 'sec6_4a.tbl1', 'h':Fdir + 'sec6_4a.tbl2', 'k':Fdir + 'sec6_4a.tbl3',
              'kp':Fdir + 'kepler_response_hires1'}



    return outflux


def run_LC(file):
    '''
    Hand a Kepler light curve (maybe with some pre-processing done)

    Call the needed methods to output ugriz light curves
    '''

    Tstar = 5000.
    Tspot = 4500.

    ffactor = measure_ampl(flux, Tstart=Tstar, dTspot=Tstar-Tspot)

    u,g,r,i,z = transform(ffactor, Tstart=Tstar, dTspot=Tstar-Tspot)

    return


# how to run code from terminal, like this:
# $ python transform.py filename.dat
if __name__ == "__main__":
    import sys
    run_LC(sys.argv[1])
