# 2MASS Reference: http://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html
#
# The 2MASS J, H and Ks relative spectral response curves (RSRs), peak-normalized to unity, derived by 
# Cohen et al. (2003). As stated by these authors, these curves "are designed to be integrated directly 
# over stellar spectra in F_lambda form, in order to calculate synthetic photometric magnitudes. The QE-based 
# component was converted to yield photon-counting RSRs by multiplying by wavelength and renormalized, 
# as described by Bessel (2000)." 


# Johnson Reference: https://en.wikipedia.org/wiki/Photometric_system (!), zero points from DCB



def getMag_hires(band,wavelength,flux,unit):
    
    import numpy as np
    import scipy.interpolate as interp
    
    
    
    '''if band=='Ks':
        bandwav,bandpass=np.loadtxt("../DATA/Ks_2MASS.txt",unpack=True)    # micron
        center=2.159
        F0= 4.283E-14        # W cm^-2 micron^-1

    elif band=='H':
        bandwav,bandpass=np.loadtxt("../DATA/H_2MASS.txt",unpack=True)
        center=1.662
        F0= 1.133E-13
    
    elif band=='J':
        bandwav,bandpass=numpy.loadtxt("../DATA/J_2MASS.txt",unpack=True)
        center=1.235
        F0= 3.129E-13'''
    
    #Units of Ks, H and J filters are microns, all other filters in Angstrom's
    
    
    if band=='Ks':
        bandwav,bandpass=np.loadtxt('filters/k_filter.txt',unpack=True)
        center=2.159        #micron
        F0=4.283E-14        #W cm^-2 micron^-1

    elif band=='H':
        bandwav,bandpass=np.loadtxt('filters/h_filter.txt',unpack=True)
        center=1.662
        F0=1.133E-13

    elif band=='J':
        bandwav,bandpass=np.loadtxt('filters/j_filter.txt',unpack=True)
        center=1.235
        F0=3.129E-13

    elif band=='U':
        bandwav,bandpass=np.loadtxt('filters/bessell_U.dat',unpack=True)
        center=0.365
        F0=4.19E-12

    elif band=='B':
        bandwav,bandpass=np.loadtxt('filters/bessell_B.dat',unpack=True)
        center=0.445
        F0=6.60E-12

    elif band=='V':
        bandwav,bandpass=np.loadtxt('filters/bessell_V.dat',unpack=True)
        center=0.551
        F0=3.61E-12

    elif band=='R':
        bandwav,bandpass=np.loadtxt('filters/bessell_R.dat',unpack=True)
        center=0.658
        F0=2.25E-12

    elif band=='I':
        bandwav,bandpass=np.loadtxt('filters/bessell_I.dat',unpack=True)
        center=0.806
        F0=1.22E-12
                                      #SDSS centers/widths/F0 from http://www.astronomy.ohio-state.edu/~martini/usefuldata.html
    elif band=='u':                   #NOTE: ugriz filters are on the AB magnitude system, KsHJUBVRI are on the Vega system
        bandwav,bandpass=np.loadtxt('filters/u_filter.txt',unpack=True)
        center=0.356         
        F0=0.8595E-11        
        
    elif band=='g':
        bandwav,bandpass=np.loadtxt('filters/g_filter.txt',unpack=True)
        center=0.483
        F0=0.4669E-11
        
    elif band=='r':
        bandwav,bandpass=np.loadtxt('filters/r_filter.txt',unpack=True)
        center=0.626
        F0=0.2780E-11
    
    elif band=='i':
        bandwav,bandpass=np.loadtxt('filters/i_filter.txt',unpack=True)
        center=0.767
        F0=0.1852E-11
    
    elif band=='z':
        bandwav,bandpass=np.loadtxt('filters/u_filter.txt',unpack=True)
        center=0.910
        F0=0.1315E-11
    
    '''elif band=='Kep':
        bandwav,bandpass=np.loadtxt()
        center=6400/10**4
        F0= '''


    filterband=np.zeros(wavelength.size)
    if band=='U'or band=='B'or band=='V'or band=='R'or band=='I'or band=='u'or band=='g'or band=='r'or band=='i'or band=='z':
        bandwav=[i/1.0E4 for i in bandwav]                     #Convert bandwav from Angstroms to microns
        
    bandinterp=interp.interp1d(bandwav,bandpass)
            #1D function between bandwav(x) and bandpass(y), y=f(x)
        
    inband=np.logical_and(wavelength>bandwav[0],wavelength<bandwav[-1] )
            #Array of 'True's where wavelength is in range of bandwav, 'False's outside of that range. Basically just an index of bandwav in wavelength array
        
    filterband[inband]=bandinterp(wavelength[inband])
            #Ignore all values of bandpass outside of filter range. Now filterband is the same shape as wavelength and flux
        
    
    if unit=='angstrom':  # convert to microns
        wavelength=wavelength/1.0E4
        flux=flux*1.0E4
    
    dwav=np.zeros(wavelength.size)
        
    dwav[0:-1]=wavelength[1:]-wavelength[0:-1]
    

    
    mag=2.5*np.log10((np.sum(F0*dwav*filterband))/(np.sum(flux*dwav*filterband)))
    
    if unit=='angstrom':
        return band,center*1.0E4,mag
    else:
        return band,center,mag 
    

