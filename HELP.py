import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
import pandas as pd

import getMag_hires
getMag_hires = getMag_hires.getMag_hires

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

def Dotter_rad( PhTemp,flux_Phoenix ):
    '''
    Reads in data from Dotter file to find stellar temperature vs radius relation used for converting flux spectra onto an absolute scale for later magnitude calculation.
    Paramters
    ---------
    PhTemp: str
        Photospheric temperature of star input by user
    
    flux_Phoenix: array
        Original spectra read in from Phoenix file, units of W/cm^2/micron
    
    Return
    ------
    Flux: array
        Converted spectra now on an absolute scale, i.e. flux recieved from 10PC, units W/cm^2/micron
    '''
    DotterData = open('models/Dotter_2014_isochrones','r')            #Opening the Dotter Data file, creating a readable table
    names = DotterData.readline().split
    DotterData.close()
    DotterData = pd.read_table('models/Dotter_2014_isochrones', sep=r'\s+', engine='python')
    
    DotterData['logAge']                                              #Look at the list of logAge data
    Age_idx = np.where(DotterData['logAge']==8.1)                     #Create an index for the data where logAge = 8.1
    TempData = 10**DotterData['logTeff'].values[Age_idx]              #Create array of temperature data, on absolute scale
    RadData = 69.63*(10**9)*10**DotterData['logR_Ro'].values[Age_idx] #Create an array of radius data, absolute scale, no ratio
    Data = list(zip(TempData,RadData))                                #Combine the temp and rad data into a tuple array so temp and rad don't become miss-ordered
    
    def getKey(item):
        return item[0]
    SortedData = sorted(Data, key=getKey)                         #Sort the tuple array by ascending temperature

    del SortedData[12]                                            #Remove the red giant data point

    Temperaturesx = np.arange(2300,12100,100).tolist()            #Creating final lists for calculations
    StrTemperaturesx = [str(i) for i in Temperaturesx]

    Radiiy = np.interp(Temperaturesx,[i[0] for i in SortedData], [i[1] for i in SortedData])
    
    PhTemp_idx = StrTemperaturesx.index(PhTemp)
    StarRad = Radiiy[PhTemp_idx]                         #StarRad is the stellar radius used for finding absolute flux

    Flux = []
    for i in range(len(flux_Phoenix)):
        y = flux_Phoenix[i]*(StarRad**2)/((3.0857*10**19)**2) #3.0857*10**19 is 10 parsecs in cm
        Flux.append(y)

    Flux = np.asarray(Flux)
    return Flux

def MagConvert( syn_phot ):
    '''
    Converts UBVRUJHKs (not ugriz) magnitude values onto the AB scale.
    Parameters
    ----------    
    syn_phot: triple-list
        Array filled with bandnames (str), centerpoints (float), and magnitudes (float)  
    
    Returns
    -------
    Filters, Centers, Magntiudes: lists
        Filled with band info, now all on AB magnitude scale
    '''
    ABtoVegaMag = [0.79,-0.09,0.02,0.21,0.45,0.91,1.39,1.85,0.91,-0.08,-0.16,-0.37,0.54]  #Conversion factors (m_AB-m_Vega)
    Filters = [i[0] for i in syn_phot]               #Breaking up syn_phot into respective lists
    Centers = [i[1] for i in syn_phot]
    Magnitudes = [i[2] for i in syn_phot]
    
    idx=0    
    for i in range(len(Filters)):
        if Filters[i]=='U' or Filters[i]=='B' or Filters[i]=='V' or Filters[i]=='R' or Filters[i]=='I' or Filters[i]=='J' or Filters[i]=='H' or Filters[i]=='Ks':
            idx=idx+1
    for i in range(idx):                             #Coverting magnitude values to AB zero-point
        Magnitudes[i] = Magnitudes[i]+ABtoVegaMag[i]
    
    syn_phot = list(zip(Filters,Centers,Magnitudes)) #Reforming the syn_phot triple
    
    def getKey(item):                                #Sorting values by ascending center-point wavelength
        return item[1]
    syn_phot = sorted(syn_phot, key=getKey)
    
    Filters = [i[0] for i in syn_phot]               #Breaking up SORTED syn_phot into respective lists
    Centers = [i[1] for i in syn_phot]
    Magnitudes = [i[2] for i in syn_phot]
    
    return Filters, Centers, Magnitudes

def V_VKPlotData(PhMax,PhMin,SpTempFrac,FillFactor):
    '''
    Routine to calculate V and (V-K) magnitudes for an input temperature range, spot temperature, and filling-factor. Primary data sets relate to temperature range as to model changes in magnitudes as a function of differing spot temperatures or filling factors (Calculated by multiple function calls).
    Parameters
    ----------
    PhMax: int
        Maximum photospheric temperature to test, must be an increment of 100 Kelvin from 2300 to 12,000
        
    PhMin: int
        Minimum photospheric temperature to test, must be an increment of 100 Kelvin from 2300 to 12,000
        
    SpTempFrac: float/array of floats
        Fractional spot temperature relative to photospheric temp
        
    FillFactor: float/array of floats
        Fractional percentage of stellar surface covered by spot
        
    Returns
    -------
    Vmags: list of floats
        V magnitudes for input parameters
        
    VKmags: list of floats
        (V-K) colors for input parameters
    '''
    PhRange = np.arange(PhMax,PhMin,-200)                              #Array of photospheric temperatures to run through
    SpRange = [str(int(round(SpTempFrac*i/100)*100)) for i in PhRange] #List of strings of spot temps as a fraction of PhTemp
    PhRange = [str(i) for i in PhRange]                                #Convert PhRange to list of strings

    FluxRange = []
    for i in range(len(PhRange)):
        wavelength,flux_Phoenix,PhTemp,SpTemp,FillFactor = make_spotmodel(PhRange[i],SpRange[i],FillFactor)
        Flux = Dotter_rad( PhTemp, flux_Phoenix )
        FluxRange.append(Flux)                 #FluxRange filled with spectra for photospheric temperatures in specified range

    Vmags = []
    for i in range(len(FluxRange)):
        band,center,mag = getMag_hires('V',wavelength,FluxRange[i],'microns')
        Vmags.append(mag)                      #Vmags filled with V magnitudes for photospheric temperatures in specified range

    Kmags = []
    for i in range(len(FluxRange)):
        band,center,mag = getMag_hires('Ks',wavelength,FluxRange[i],'microns')
        Kmags.append(mag)                     #Kmags filled with K magnitudes for photospheric temperatures in specified range

    VKmags = [i-j for i,j in zip(Vmags,Kmags)]#VKmags filled with (V-K) colors for photospheric temperatures in specified range
    
    return(Vmags,VKmags)

def DeltaVfunc(x,y):
    '''
    General function to calculate the difference between the V vs (V-K) Kamaii isochrone and V magnitudes of data set.
    Parameters
    ----------
    x: list of floats
        (V-K) color values of data set
    y: list of floats
        V magnitudes of data set
    Returns
    -------
    DeltaV_VK: list of two lists
        (V-K) magnitudes stored in DeltaV_VK[0]
        DeltaV magnitudes stored in DeltaV_VK[1]
    '''
    KamaiIsochrone = open('data/Kamai_Isochrone.tbl','r')
    names = KamaiIsochrone.readline().split
    KamaiIsochrone.close()
    KamaiIsochrone = pd.read_table('data/Kamai_Isochrone.tbl',sep=None,engine='python')
    
    KamaiV = [i-5.67 for i in KamaiIsochrone['V']]
    KamaiVK = KamaiIsochrone['V - K']
    
    VFUNC = np.interp(x,KamaiVK,KamaiV)
    
    DeltV = []                    #Vertical separation between data V magnitudes and Kamai isochrone
    for i in range(len(y)):
        z = y[i] - VFUNC[i]
        DeltV.append(z)
    
    DeltaV_VK = list(zip(x,DeltV))
    
    return DeltaV_VK

def PleiadesDeltaV_VK(x,y):
    '''
    Returns Delta V and (V-K) arrays for Pleiades data with respect to Kamai isochrone. Could be modified to a more general function that inputs more filters, currently limited up to V,B,I,K due to Kamaii data.
    Parameters
    ----------
    x: str
        Currently of no use, would use to specify filters
    y: str
        ""
    
    Returns
    -------
    PleiaVK: list of floats
        (V-K) colors, used for x axis of CMD    
    PleiaDeltaV: list of floats
        Delta V differences SPECIFICALLY between Pleiades V magnitudes and Kamai isochrone
    '''
    
    PleiadesData = open('data/PleiadesColors.tbl','r')
    names = PleiadesData.readline().split
    PleiadesData.close()
    PleiadesData = pd.read_table('data/PleiadesColors.tbl',sep=r'\s+',engine='python')

    PleiaV = PleiadesData['V_mag']
    PleiaV = [i-5.67 for i in PleiaV] #Puts magnitudes on absolute scale with m-M=5.67 distance modulus
    
    PleiaK = PleiadesData['K_mag']
    PleiaK = [i-5.67 for i in PleiaK]
    
    PleiaVK = [i-j for i,j in zip(PleiaV,PleiaK)]
    
    PleiaVVK = list(zip(PleiaVK,PleiaV))
    PleiaVVK = [[x,y] for x,y in PleiaVVK if x>=1 and x<=6.4]
    
    PleiaVK = [i[0] for i in PleiaVVK]
    PleiaV = [i[1] for i in PleiaVVK]
    
    DeltaV_VK = DeltaVfunc(PleiaVK,PleiaV)
    
    DeltaV_VK = [[x,y] for x,y in DeltaV_VK if y>-0.375]

    PleiaVK = [i[0] for i in DeltaV_VK]
    PleiaDeltaV = [i[1] for i in DeltaV_VK]
    
    return PleiaVK,PleiaDeltaV

def CMDData(x,y,z=None):
    '''
    Returns magnitude data for specified filter bands from Pleiades data table.
    Parameters
    ----------
    x: str
        x-axis for CMD
    y: str
        y-axis for CMD
    z: str
        color for CMD if wanted
    Returns
    -------
    Pleiax: float list
        List of x axis magnitudes
    Pleiay: float list
        List of y axis magnitudes
    '''
    PleiadesData = open('data/PleiadesColors.tbl','r')
    names = PleiadesData.readline().split
    PleiadesData.close()
    PleiadesData = pd.read_table('data/PleiadesColors.tbl',sep=r'\s+',engine='python')
    
    Pleiax=PleiadesData[x+'_mag']
    Pleiax = [i-5.67 for i in Pleiax]
    Pleiay = PleiadesData[y+'_mag']
    Pleiay = [i-5.67 for i in Pleiay]
    
    return Pleiax,Pleiay