'''
Code to do transformation from offset in 1 filter system to offset in another

Doing with mindset (currently) of:
given starspot modulation from Kepler, a spot delta T, and a stellar T_eff
 predict the ugrizJHK response
'''


import numpy as np
import matplotlib.pyplot as plt
from astropy.analytic_functions import blackbody_lambda


def transform(dflux, Fin='kp', Fout='g', Tstar=5777.0, dTspot=500.):
    '''

    :param dmag:
    :param Fin:
    :param Fout:
    :return:
    '''

    Fdir = 'filters/'
    Ffiles = {'u':Fdir + 'u.dat.txt', 'g':Fdir + 'g.dat.txt', 'r':Fdir + 'r.dat.txt',
              'i':Fdir + 'i.dat.txt', 'z':Fdir + 'z.dat.txt',
              'j':Fdir + 'sec6_4a.tbl1', 'h':Fdir + 'sec6_4a.tbl2', 'k':Fdir + 'sec6_4a.tbl3',
              'kp':Fdir + 'kepler_response_hires1'}


    # make BB w/ star Temp

    # make BB w/ spot Temp

    # take input band's flux, convolve output flux to proper band

    return outflux