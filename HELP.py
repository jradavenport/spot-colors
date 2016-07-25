import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits


def get_model( temp ):
    '''
    Function opens Phoenix file of specified temperature in 100 Kelvin increments from 2300 to 12000 Kelvin.
    Parameters
    ----------
    temp: str
        Desired temperature in 100K increments from 2300 to 12,000 Kelvin
        
    Returns
    -------
    wavelength: array
        Associated wavelengths for flux measurements from Phoenix FITS file, in Angstroms
        
    flux: array
        Flux measurements from Phoenix FITS file, in erg/s/cm^2
        
    temp: str
        See Parameters
    '''
    model_file = 'PHOENIX-ACES-AGSS-COND-2011_R10000FITS_Z-0.0/lte0' + temp + '-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
    hdulist = fits.open(model_file)
    flux = hdulist[0].data
    wavelength = np.exp(hdulist[0].header[('CRVAL1')]+hdulist[0].header[('CDELT1')]*np.arange(0,212027))/(10**4)
    return wavelength, flux, hdulist
    
def make_spotmodel( PhTemp, SpTemp, FillFactor ):
    '''
    Funtion models a spotted star through the get_model function by adding two different temperatured spectra
    and combining with a certain filling factor.
    Parameters
    ----------
    PhTemp: str
        Desired photospheric temperature of the star, in 100 Kelvin increments from 2300 to 12,000
    
    SpTemp: str
        Desired temperature of spot, same scale as PhTemp
    
    FillFactor: int
        Fraction of the star's surface covered by spot, decimal from 0 to 1
        
    Returns
    -------
    wavelength: array
        Corresponding wavelengths of spectral data, units of microns
    
    tot_flux: array
        Combined fluxes of the photospheric and spot spectra, units of W/cm^2/micron
    
    PhTemp, SpTemp, FillFactor: str
        Input values
    '''
    PhWavelength, PhFlux, hdulist1=get_model( PhTemp )
    SpWavelength, SpFlux, hdulist2=get_model( SpTemp )
    flux_Phoenix = (10**-11)*((1-FillFactor)*PhFlux + FillFactor*SpFlux)
    wavelength = (np.exp(hdulist1[0].header[('CRVAL1')]+hdulist1[0].header[('CDELT1')]*np.arange(0,212027)))/(10**4)
    return wavelength, flux_Phoenix, PhTemp, SpTemp, FillFactor