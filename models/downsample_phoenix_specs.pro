FUNCTION PLANCK, temp, wave
COMMON constants, kerg, kev, herg, me, mp, ang, ev, ephot, c, Msun, grav, Lsun,Rsun, parsec, AU, solarmagbc, stefboltz, thompson

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                         ;
;       This function computes the value of the Planck function in cgs    ;
;       units given a temperature (K) and wavelength (cm) of interest.    ;
;       The value of the function is returned as a floating point value.  ;
;                                                                         ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
kerg = 1.38065D-16
kev  = 8.61734D-5
herg = 6.62607D-27
me   = 9.10939D-28
mp   = 1.67262D-24
ang  = 1D-8
ev   = 1.60219D-12
ephot= 12398.55
c    = 2.99792D10
Msun = 1.989D33
grav = 6.67259D-8
Lsun = 3.826D33
Rsun = 6.9599D10
parsec = 3.0857D18
AU   = 1.496D13
solarmagbc = 4.65
stefboltz = 5.6705D-5
thompson = 6.648D-25


;calculate the planck function for the temperature and wavelength supplied 
blambda = 2*herg*c*c/((wave^5)*(EXP(herg*c/(wave*kerg*temp))-1))
 
RETURN, blambda
 
END

;---------------------------

FUNCTION IReqw, lambda, flux, linecent, linewidth, cont1, cont1width, cont2, cont2width, PLOT = PLOT

;interpolate the continuum over the line
incont1 = WHERE(lambda GT cont1-cont1width/2. AND lambda LT cont1+cont1width/2., n_incont1)
incont2 = WHERE(lambda GT cont2-cont2width/2. AND lambda LT cont2+cont2width/2., n_incont2)
inline = WHERE(lambda GT linecent-linewidth/2. AND lambda LT linecent+linewidth/2., n_inline)

IF n_incont1 GT 1 AND n_incont2 GT 1 AND n_inline GT 1 THEN BEGIN

    ;find the median flux in the continuum regions
    mediancont1 = MEDIAN(flux[incont1])
    mediancont2 = MEDIAN(flux[incont2])
    
    ;interpolate between continuum
    ;regions through the line
    continuum = INTERPOL([mediancont1,mediancont2],[cont1,cont2],lambda[inline])

    ;find the flux absorbed in the line
    absorbed_flux = TSUM(lambda[inline], continuum - flux[inline])

    ;find the mean flux in the line
    line_flux = TSUM(lambda[inline], flux[inline])/(lambda[inline[n_inline-1]]-lambda[inline[0]])

    ;divide by the continuum flux density in the center of the line
    distance_to_linecenter = linecent-lambda[inline]
    closest = MIN(ABS(distance_to_linecenter), closest_index)
    eqw = absorbed_flux*1./continuum[closest_index]
    flux_ratio = line_flux*1./continuum[closest_index]

    ;make the index a magnitude
    index = -2.5*ALOG10(flux_ratio)

    IF KEYWORD_SET(PLOT) THEN BEGIN
        plotstart = cont1-4*cont1width
        plotend = cont2+4*cont2width

        onplot = WHERE(lambda GT plotstart AND lambda LT plotend)
        hiplot = MAX(flux[onplot])
        loplot = MIN(flux[onplot])

        PLOT, lambda, flux, XRANGE = [plotstart,plotend],YRANGE = [0.8*loplot,1.2*hiplot], PSYM = 10, /XSTY, /YSTY
        OPLOT, lambda[inline], continuum, LINESTYLE = 2
        OPLOT, [linecent+linewidth/2.,linecent+linewidth/2.], [-1000.,1000.], LINESTYLE = 3
        OPLOT, [linecent-linewidth/2.,linecent-linewidth/2.], [-1000.,1000.], LINESTYLE = 3
        OPLOT, [cont1-cont1width/2.,cont1-cont1width/2.], [-1000.,1000.], LINESTYLE = 4
        OPLOT, [cont1+cont1width/2.,cont1+cont1width/2.], [-1000.,1000.], LINESTYLE = 4
        OPLOT, [cont2-cont2width/2.,cont2-cont2width/2.], [-1000.,1000.], LINESTYLE = 4
        OPLOT, [cont2+cont2width/2.,cont2+cont2width/2.], [-1000.,1000.], LINESTYLE = 4
        OPLOT, [lambda[inline[closest_index]],lambda[inline[closest_index]]], [continuum[closest_index], continuum[closest_index]], PSYM = 4
        
        ;wait = GET_KBRD(1)

    ENDIF

ENDIF ELSE BEGIN
    eqw = -99.
    index = -99.
ENDELSE

RETURN, eqw*10000. ;index ;

END
;-----------------------------------------------
FUNCTION IRratio, lambda, flux, numcent, numwidth, denomcent, denomwidth, PLOT = PLOT

innum = WHERE(lambda GT numcent-numwidth/2. AND lambda LT numcent+numwidth/2., n_num)
indenom = WHERE(lambda GT denomcent-denomwidth/2. AND lambda LT denomcent+denomwidth/2., n_denom)

IF n_num GT 1 AND n_denom GT 1 THEN BEGIN

    ratio = MEAN(flux[innum]/flux[indenom])

    IF KEYWORD_SET(PLOT) THEN BEGIN
        plotrange = MAX([lambda[innum],lambda[indenom]]) - MIN([lambda[innum],lambda[indenom]])
        plotstart = MIN([lambda[innum],lambda[indenom]]) - plotrange/2.
        plotend = MAX([lambda[innum],lambda[indenom]]) + plotrange/2.

        onplot = WHERE(lambda GT plotstart AND lambda LT plotend)
        hiplot = MAX(flux[onplot])
        loplot = MIN(flux[onplot])

        PLOT, lambda, flux, XRANGE = [plotstart,plotend],YRANGE = [0.8*loplot,1.2*hiplot], PSYM = 10, /XSTY, /YSTY
        OPLOT, [numcent-numwidth/2.,numcent-numwidth/2.], [-1000.,1000.], LINESTYLE = 4
        OPLOT, [numcent+numwidth/2.,numcent+numwidth/2.], [-1000.,1000.], LINESTYLE = 4
        OPLOT, [denomcent-denomwidth/2.,denomcent-denomwidth/2.], [-1000.,1000.], LINESTYLE = 4
        OPLOT, [denomcent+denomwidth/2.,denomcent+denomwidth/2.], [-1000.,1000.], LINESTYLE = 4
        
        wait = GET_KBRD(1)

    ENDIF

ENDIF ELSE BEGIN

    ratio = -99.

ENDELSE

RETURN, ratio

END

;-----------------------------------------------
FUNCTION measureIRindices, lambda, flux
 
;define an array where I can store the spectral indices
indices = FLTARR(30)

;start measuring indices

;measure Mg I @ 1.48 microns in the H band
indices[0] = IREQW(lambda, flux, 1.4892, 0.006, 1.482, 0.003, 1.494, 0.003)

;measure Mg I @ 1.50 microns in the H band
indices[1] = IREQW(lambda, flux, 1.504, 0.004, 1.498, 0.0045, 1.5095, 0.0045)

;measure K I @ 1.52 microns in the H band
indices[2] = IREQW(lambda, flux, 1.5172, 0.004, 1.5105, 0.004, 1.523, 0.004)

;measure Mg I @ 1.576 
indices[3] = IREQW(lambda, flux, 1.576, 0.004, 1.57, 0.004, 1.58, 0.004)

;measure a simple ratio using two equivalent widths with the same
;continuum regions...
indices[4] = -IREQW(lambda, flux, 1.575, 0.002, 1.57, 0.004, 1.5815, 0.003)+IREQW(lambda, flux, 1.5775, 0.003, 1.57, 0.004, 1.5815, 0.003)
;indices[4] = IRRATIO(lambda, flux, 1.575, 0.002, 1.5775, 0.003, /plot)

;measure Si I @ 1.59 microns in the H band
indices[5] = IREQW(lambda, flux, 1.59, 0.005, 1.586, 0.003, 1.594, 0.003)

;measure CO @ 1.62 microns in the H band
indices[6] = IREQW(lambda, flux, 1.61975, 0.0045, 1.615, 0.003, 1.628, 0.003)

;measure CO @ 1.66 (inspired by Meyer et al. 1998)
indices[7] = IREQW(lambda, flux, 1.6625, 0.004, 1.6565, 0.003, 1.67, 0.003)

;measure Al @ 1.67 microns in the H band
indices[8] = IREQW(lambda, flux, 1.674, 0.006, 1.659, 0.003, 1.678, 0.002)

;measure HI @ 1.68 microns in the H band
indices[9] = IREQW(lambda, flux, 1.681, 0.004, 1.6775, 0.003, 1.6855, 0.003)

;measure OH in H band (inspired by Meyer et al. 1998)
indices[10] = IREQW(lambda, flux, 1.686, 0.0065, 1.6785, 0.003, 1.695, 0.003)

;Measure bump next to 1.71 Mg line
indices[11] = IREQW(lambda, flux, 1.7075, 0.003, 1.704, 0.003, 1.7145, 0.003)

;Measure Mg I @ 1.71 microns in the H band
indices[12] = IREQW(lambda, flux, 1.7115, 0.003, 1.704, 0.003, 1.7145, 0.003)

;measure a continuum slope in the H band
indices[13] = IRRATIO(lambda, flux, 1.605, 0.02, 1.69, 0.02)

;measure water absorption in the H band
indices[14] = IRRATIO(lambda, flux, 1.69, 0.02, 1.77, 0.02)

;Measure Ca I @ 1.982 microns in the H band
indices[15] = IREQW(lambda, flux, 1.982, 0.013, 1.9676, 0.005, 1.99775, 0.005)

;measure Mg I @ 2.1066 + Al I 2.1099 microns in the K band 
indices[16] = IREQW( lambda, flux, 2.1082, 0.006, 2.102, 0.004, 2.113, 0.004)

;measure Mg I @ 2.1066 microns in the K band 
indices[17] = IREQW( lambda, flux, 2.1065, 0.003, 2.102, 0.004, 2.113, 0.004)

;measure Al I @ 2.1099 microns in the K band 
indices[18] = IREQW( lambda, flux, 2.1097, 0.003, 2.102, 0.004, 2.113, 0.004)

;measure Al I @ 2.1168 microns in the K band
indices[19] = IREQW( lambda, flux, 2.1168, 0.003, 2.1145, 0.002, 2.1195, 0.002)

;measure Brackett Gamma @ 2.165 microns in the K band (using
;bandpasses from Kleinmann & Hall 86, but with a modified red continuum)
indices[20] = IREQW( lambda, flux, 2.16625, 0.0047, 2.155, 0.0044, 2.175, 0.0042)

;measure Na I doublet @ 2.21 microns in the K band (using bandpasses
;from Rojas-Ayala et al. 2012)
indices[21] = IREQW( lambda, flux, 2.207, 0.005, 2.1965, 0.004, 2.2175, 0.004)

;measure asymmetry in NaI profile for giant/dwarf separation
indices[22] = IREQW( lambda, flux, 2.206, 0.0035, 2.195, 0.004, 2.217, 0.006)
indices[23] = IREQW( lambda, flux, 2.2095, 0.0035, 2.195, 0.004, 2.217, 0.006);/IREQW( lambda, flux, 2.206, 0.0035, 2.195, 0.004, 2.217, 0.006)

;measure Ca I @ 2.26 microns in the K band (using bandpasses
;from Royas-Ayala et al. 2012)
indices[24] = IREQW( lambda, flux, 2.2635, 0.0055, 2.2545, 0.0045, 2.273, 0.0045)

;measure Mg I @ 2.28 microns in the K band (using bandpasses
;from Ivanov et al. 04)
indices[25] = IREQW( lambda, flux, 2.2806, 0.0052, 2.275, 0.006, 2.2857, 0.005)

;measure CO @ 2.3 microns in the K band (using bandpasses 
;from Ivanov et al. 06, plus a long wavelength anchor I invented)
indices[26] = IREQW(lambda, flux, 2.30375, 0.0225, 2.288, 0.007, 2.3185, 0.004)

;measure a continuum slope in the K band for constructing H20_K2
indices[27] = IRRATIO(lambda, flux, 2.08, 0.02, 2.245, 0.02)

;measure water absorption in the K band for constructing H20_K2
indices[28] = IRRATIO(lambda, flux, 2.245, 0.02, 2.37, 0.02)

;measure co feature at 2.345 microns
indices[29] = IREQW(lambda, flux, 2.3455, 0.003, 2.3425, 0.003, 2.349, 0.003)

RETURN, indices

END

;-------
;----------------------------------------------------------------------

function IRresponse,ilambda, ivac_pass, lambda
;this interpolates over the 2MASS filter curves and matches the number
;of elements in the data

response=INTERPOL(ivac_pass, ilambda, lambda)

return,response

end

;----------------------------------------------------------------------

FUNCTION IRfilter, lambda, flux, jlambda, jcurve, hlambda, hcurve, klambda, kcurve, W1lambda, W1curve, W2lambda, W2curve, W3lambda, W3curve, W4lambda, W4curve

;PRINT, lambda
;PRINT, W1lambda

;these indices define where the data overlap the curves. 
jind = WHERE(lambda GE MIN(jlambda) AND lambda LE MAX(jlambda))
hind = WHERE(lambda GE MIN(hlambda) AND lambda LE MAX(hlambda))
kind = WHERE(lambda GE MIN(klambda) AND lambda LE MAX(klambda))
;w1ind = WHERE(lambda GE MIN(W1lambda) AND lambda LE MAX(W1lambda))
;w2ind = WHERE(lambda GE MIN(W2lambda) AND lambda LE MAX(W2lambda))
;w3ind = WHERE(lambda GE MIN(W3lambda) AND lambda LE MAX(W3lambda))
;w4ind = WHERE(lambda GE MIN(W4lambda) AND lambda LE MAX(W4lambda))

;PRINT, w1ind

;;these are testing where the filters overlap the data (comment out in
;;the ideal case of having full wavelength coverage)
;jresind = WHERE(jlambda GE MIN(lambda) AND jlambda LE MAX(lambda))
;hresind = WHERE(hlambda GE MIN(lambda) AND hlambda LE MAX(lambda))
;kresind = WHERE(klambda GE MIN(lambda) AND klambda LE MAX(lambda))

;do the photometry for each band individually, so we can test for
;bands with missing data.  Only comment the J band to avoid redundancy

    ;compute the response in this band
    jres = IRRESPONSE(jlambda, jcurve, lambda[jind])
    hres = IRRESPONSE(hlambda, hcurve, lambda[hind])
    kres = IRRESPONSE(klambda, kcurve, lambda[kind])
;    w1res = IRRESPONSE(w1lambda, w1curve, lambda[w1ind])
;    w2res = IRRESPONSE(w2lambda, w2curve, lambda[w2ind])
;    w3res = IRRESPONSE(w3lambda, w3curve, lambda[w3ind])
;    w4res = IRRESPONSE(w4lambda, w4curve, lambda[w4ind])

    ;take a flux at each wavelength in units of ergs/cm^2-sec-Ang,
    ;convolve with filter response, use the fact that dn=dl^2/c where dl
    ;is a change in wavelength and dn is a change in frequency to convert
    ;the 1/A to 1/Hz (remember that 2.998d18 is c in angstroms per second
    ;to make the units work out right), then use the fact that 1 Jy =
    ;10^-23 ergs/sec/cm^2/Hz to get the flux at each wavelength in janskys.
    jflux = DOUBLE(jres*flux[jind]*(lambda[jind]^2/2.998d18)*1.d23)
    hflux = DOUBLE(hres*flux[hind]*(lambda[hind]^2/2.998d18)*1.d23)
    kflux = DOUBLE(kres*flux[kind]*(lambda[kind]^2/2.998d18)*1.d23)
;    w1flux = DOUBLE(w1res*flux[w1ind]*(lambda[w1ind]^2/2.998d18)*1.d23)
;    w2flux = DOUBLE(w2res*flux[w2ind]*(lambda[w2ind]^2/2.998d18)*1.d23)
;    w3flux = DOUBLE(w3res*flux[w3ind]*(lambda[w3ind]^2/2.998d18)*1.d23)
;    w4flux = DOUBLE(w4res*flux[w4ind]*(lambda[w4ind]^2/2.998d18)*1.d23)

    ;integrate using trapezoidal rule to get the total flux (change 
    ;the ordinate into frequency space by using c [in angstroms per 
    ;second] / wavelength to get Hz) and then divide by the "effective" 
    ;filter width (which is essentially the width a perfect [response = 1]
    ;boxcar filter would have to have to get the same amount of flux
    ;through as our real filter did) to get the flux density in Janskys 
    jjansk = ABS(TSUM(2.998d18/lambda[jind],jflux)/(2.998d18*TSUM(lambda[jind],jres)/12350.^2))
    hjansk = ABS(TSUM(2.998d18/lambda[hind],hflux)/(2.998d18*TSUM(lambda[hind],hres)/16620.^2))
    kjansk = ABS(TSUM(2.998d18/lambda[kind],kflux)/(2.998d18*TSUM(lambda[kind],kres)/21590.^2)) ;6205
;    w1jansk = ABS(TSUM(2.998d18/lambda[w1ind],w1flux)/(2.998d18*TSUM(lambda[w1ind],w1res)/33526.^2)) ;6205
;    w2jansk = ABS(TSUM(2.998d18/lambda[w1ind],w1flux)/(2.998d18*TSUM(lambda[w1ind],w1res)/46028.^2)) ;6205
;    w3jansk = ABS(TSUM(2.998d18/lambda[w1ind],w1flux)/(2.998d18*TSUM(lambda[w1ind],w1res)/115608.^2)) ;6205
;    w4jansk = ABS(TSUM(2.998d18/lambda[w1ind],w1flux)/(2.998d18*TSUM(lambda[w1ind],w1res)/220883.^2)) ;6205

fluxes = [jjansk, hjansk, kjansk] ;, w1jansk, w2jansk, w3jansk, w4jansk]

RETURN,fluxes

END

;-------------------------------------------------------------------

function jansk2IRmag, janskys

;calculate using 2MASS zeropoints (in janskys) as reported by
;Cohen et al 2003

Jmag = -2.5*ALOG10(janskys[0]/1594.)
Hmag = -2.5*ALOG10(janskys[1]/1024.)
Kmag = -2.5*ALOG10(janskys[2]/666.7)
;W1mag = -2.5*ALOG10(janskys[3]/307.)
;W2mag = -2.5*ALOG10(janskys[4]/171.)
;W3mag = -2.5*ALOG10(janskys[5]/30.3)
;W4mag = -2.5*ALOG10(janskys[6]/8.32)

mags=[Jmag, Hmag, Kmag] ;, W1mag, W2mag, W3mag, W4mag]

return, mags

end

PRO load_all_filters, ulambda, ucurve, glambda, gcurve, rlambda, rcurve, ilambda, icurve, zlambda, zcurve, $
ubesslambda, ubesscurve, bbesslambda, bbesscurve, vbesslambda, vbesscurve, rbesslambda, rbesscurve, ibesslambda, ibesscurve, $
int_wfc_u_lambda, int_wfc_u_trans, int_wfc_g_lambda, int_wfc_g_trans, int_wfc_r_lambda, int_wfc_r_trans, $
int_wfc_i_lambda, int_wfc_i_trans, int_wfc_z_lambda, int_wfc_z_trans, $
Asahi_u_lambda, Asahi_u_trans, Asahi_g_lambda, Asahi_g_trans, Asahi_r_lambda, Asahi_r_trans, Asahi_i_lambda, Asahi_i_trans, Asahi_z_lambda, Asahi_z_trans, $
Andover_u_lambda, Andover_u_trans, Andover_g_lambda, Andover_g_trans, Andover_r_lambda, Andover_r_trans, Andover_i_lambda, Andover_i_trans, Andover_z_lambda, Andover_z_trans, jlambda, jcurve, hlambda, hcurve, klambda, kcurve, W1lambda, W1curve, W2lambda, W2curve, W3lambda, W3curve, W4lambda, W4curve, Mould_r_lambda, Mould_r_trans

;read in all the filters

;the actual filters used in the kcorrect routine.  note that these
;don't distinguish between psf and extended source curves, but the
;difference between these and the psf curves above is less than their
;difference from extended sources curves, so I assume these are good
;for psf sources...
READCOL,'~/Dropbox/idl/data/SDSS/sdss_u.txt',name, ulambda,uvac,uzenith, ucurve,f='A,I,f,f,f',/SILENT
READCOL,'~/Dropbox/idl/data/SDSS/sdss_g.txt',name, glambda,uvac,uzenith, gcurve,f='A,I,f,f,f',/SILENT
READCOL,'~/Dropbox/idl/data/SDSS/sdss_r.txt',name, rlambda,uvac,uzenith, rcurve,f='A,I,f,f,f',/SILENT
READCOL,'~/Dropbox/idl/data/SDSS/sdss_i.txt',name, ilambda,uvac,uzenith, icurve,f='A,I,f,f,f',/SILENT
READCOL,'~/Dropbox/idl/data/SDSS/sdss_z.txt',name, zlambda,uvac,uzenith, zcurve,f='A,I,f,f,f',/SILENT

;read in bessell filters too
READCOL, '~/Dropbox/idl/data/bessell/bessell_U.dat', ubesslambda, ubesscurve, f='f,f', /SILENT
READCOL, '~/Dropbox/idl/data/bessell/bessell_B.dat', bbesslambda, bbesscurve, f='f,f', /SILENT
READCOL, '~/Dropbox/idl/data/bessell/bessell_V.dat', vbesslambda, vbesscurve, f='f,f', /SILENT
READCOL, '~/Dropbox/idl/data/bessell/bessell_R.dat', rbesslambda, rbesscurve, f='f,f', /SILENT
READCOL, '~/Dropbox/idl/data/bessell/bessell_I.dat', ibesslambda, ibesscurve, f='f,f', /SILENT

;read in INT_WFC filters too
READCOL, '~/Dropbox/idl/data/INT_WFC/intwfc_u.dat', int_wfc_u_lambda, int_wfc_u_trans, f='f,f', /SILENT
READCOL, '~/Dropbox/idl/data/INT_WFC/intwfc_g.dat', int_wfc_g_lambda, int_wfc_g_trans, f='f,f', /SILENT
READCOL, '~/Dropbox/idl/data/INT_WFC/intwfc_r.dat', int_wfc_r_lambda, int_wfc_r_trans, f='f,f', /SILENT
READCOL, '~/Dropbox/idl/data/INT_WFC/intwfc_i.dat', int_wfc_i_lambda, int_wfc_i_trans, f='f,f', /SILENT
READCOL, '~/Dropbox/idl/data/INT_WFC/intwfc_z.dat', int_wfc_z_lambda, int_wfc_z_trans, f='f,f', /SILENT

;read in the LMI QE curve
LMI_QE_lambdas = [3500, 4000, 5000, 6500, 9000]
LMI_QE = [42.9, 90.1, 96.8, 92.9, 59.0]

;Interp QE onto Bessell filters
uQEonbess = INTERPOL(LMI_QE, LMI_QE_lambdas, ubesslambda)
ubesscurve  = ubesscurve * uQEonbess/100.
bQEonbess = INTERPOL(LMI_QE, LMI_QE_lambdas, bbesslambda)
bbesscurve  = bbesscurve * bQEonbess/100.
vQEonbess = INTERPOL(LMI_QE, LMI_QE_lambdas, vbesslambda)
vbesscurve  = vbesscurve * vQEonbess/100.
rQEonbess = INTERPOL(LMI_QE, LMI_QE_lambdas, rbesslambda)
rbesscurve  = rbesscurve * rQEonbess/100.
iQEonbess = INTERPOL(LMI_QE, LMI_QE_lambdas, ibesslambda)
ibesscurve  = ibesscurve * iQEonbess/100.

;Interp QE onto INT_WFC filters
uQEonint = INTERPOL(LMI_QE, LMI_QE_lambdas, int_wfc_u_lambda)
int_wfc_u_trans  = int_wfc_u_trans * uQEonint/100.
gQEonint = INTERPOL(LMI_QE, LMI_QE_lambdas, int_wfc_g_lambda)
int_wfc_g_trans  = int_wfc_g_trans * gQEonint/100.
rQEonint = INTERPOL(LMI_QE, LMI_QE_lambdas, int_wfc_r_lambda)
int_wfc_r_trans  = int_wfc_r_trans * rQEonint/100.
iQEonint = INTERPOL(LMI_QE, LMI_QE_lambdas, int_wfc_i_lambda)
int_wfc_i_trans  = int_wfc_i_trans * iQEonint/100.
zQEonint = INTERPOL(LMI_QE, LMI_QE_lambdas, int_wfc_z_lambda)
int_wfc_z_trans  = int_wfc_z_trans * zQEonint/100.

;read in mould R for PTF
READCOL, '~/Dropbox/idl/data/Mould_R.dat', mould_r_lambda, mould_r_trans, f='f,f', /SILENT
mould_r_lambda = REVERSE(mould_r_lambda)*10.
mould_r_trans = REVERSE(mould_r_trans)
rQEonmould = INTERPOL(LMI_QE, LMI_QE_lambdas, mould_r_lambda)
mould_r_trans = mould_r_trans* rQEonmould

;read Asahi interference filters
READCOL, '~/Dropbox/Lowell\ Docs/DCT\ info/LMI_quotes/Asahi_Quote+Specs/u-band.txt', Asahi_u_lambda, Asahi_u_trans
READCOL, '~/Dropbox/Lowell\ Docs/DCT\ info/LMI_quotes/Asahi_Quote+Specs/g-band.txt', Asahi_g_lambda, Asahi_g_trans
READCOL, '~/Dropbox/Lowell\ Docs/DCT\ info/LMI_quotes/Asahi_Quote+Specs/r-band.txt', Asahi_r_lambda, Asahi_r_trans
READCOL, '~/Dropbox/Lowell\ Docs/DCT\ info/LMI_quotes/Asahi_Quote+Specs/i-band.txt', Asahi_i_lambda, Asahi_i_trans
READCOL, '~/Dropbox/Lowell\ Docs/DCT\ info/LMI_quotes/Asahi_Quote+Specs/z-band.txt', Asahi_z_lambda, Asahi_z_trans
;READCOL, 'Asahi_Quote+Specs/Y-band.txt', Y_lambda, Y_trans

Asahi_u_lambda = Asahi_u_lambda*10.
Asahi_g_lambda = Asahi_g_lambda*10.
Asahi_r_lambda = Asahi_r_lambda*10.
Asahi_i_lambda = Asahi_i_lambda*10.
Asahi_z_lambda = Asahi_z_lambda*10.

;Interp QE onto filters
uQEonAsahi = INTERPOL(LMI_QE, LMI_QE_lambdas, Asahi_u_lambda)
Asahi_u_trans = REVERSE(Asahi_u_trans/100. * uQEonAsahi/100.)
Asahi_u_lambda = REVERSE(Asahi_u_lambda)
gQEonAsahi = INTERPOL(LMI_QE, LMI_QE_lambdas, Asahi_g_lambda)
Asahi_g_trans = REVERSE(Asahi_g_trans/100. * gQEonAsahi/100.)
Asahi_g_lambda = REVERSE(Asahi_g_lambda)
rQEonAsahi = INTERPOL(LMI_QE, LMI_QE_lambdas, Asahi_r_lambda)
Asahi_r_trans = REVERSE(Asahi_r_trans/100. * rQEonAsahi/100.)
Asahi_r_lambda = REVERSE(Asahi_r_lambda)
iQEonAsahi = INTERPOL(LMI_QE, LMI_QE_lambdas, Asahi_i_lambda)
Asahi_i_trans = REVERSE(Asahi_i_trans/100. * iQEonAsahi/100.)
Asahi_i_lambda = REVERSE(Asahi_i_lambda)
zQEonAsahi = INTERPOL(LMI_QE, LMI_QE_lambdas, Asahi_z_lambda)
Asahi_z_trans = REVERSE(Asahi_z_trans/100. * zQEonAsahi/100.)
Asahi_z_lambda = REVERSE(Asahi_z_lambda)

;read Andover Filters
READCOL, '~/Dropbox/Lowell\ Docs/DCT\ info/LMI_quotes/AndoverQuote+Specs/Sloan_u.txt', Andover_u_lambda, Andover_u_trans
READCOL, '~/Dropbox/Lowell\ Docs/DCT\ info/LMI_quotes/AndoverQuote+Specs/Sloan_g.txt', Andover_g_lambda, Andover_g_trans
READCOL, '~/Dropbox/Lowell\ Docs/DCT\ info/LMI_quotes/AndoverQuote+Specs/Sloan_r.txt', Andover_r_lambda, Andover_r_trans;, Andover_r_reflect
READCOL, '~/Dropbox/Lowell\ Docs/DCT\ info/LMI_quotes/AndoverQuote+Specs/Sloan_i.txt', Andover_i_lambda, Andover_i_trans;, Andover_i_reflect
READCOL, '~/Dropbox/Lowell\ Docs/DCT\ info/LMI_quotes/AndoverQuote+Specs/Sloan_z.txt', Andover_z_lambda, Andover_z_trans;, Andover_z_reflect
;READCOL, 'AndoverQuote+Specs/Y-band.txt', Y_lambda, Y_trans

Andover_u_lambda = Andover_u_lambda*10.
Andover_g_lambda = Andover_g_lambda*10.
Andover_r_lambda = Andover_r_lambda*10.
Andover_i_lambda = Andover_i_lambda*10.
Andover_z_lambda = Andover_z_lambda*10.

;Interp QE onto filters
uQEonAndover = INTERPOL(LMI_QE, LMI_QE_lambdas, Andover_u_lambda)
Andover_u_trans = Andover_u_trans/100. * uQEonAndover/100.
gQEonAndover = INTERPOL(LMI_QE, LMI_QE_lambdas, Andover_g_lambda)
Andover_g_trans = Andover_g_trans/100. * gQEonAndover/100.
rQEonAndover = INTERPOL(LMI_QE, LMI_QE_lambdas, Andover_r_lambda)
Andover_r_trans = Andover_r_trans/100. * rQEonAndover/100.
iQEonAndover = INTERPOL(LMI_QE, LMI_QE_lambdas, Andover_i_lambda)
Andover_i_trans = Andover_i_trans/100. * iQEonAndover/100.
zQEonAndover = INTERPOL(LMI_QE, LMI_QE_lambdas, Andover_z_lambda)
Andover_z_trans = Andover_z_trans/100. * zQEonAndover/100.

;read in 2MASS filter curves.
READCOL,'/Users/covey/Dropbox/idl/data/twomass/twomass_J.par', FORMAT = '(X, F, F)', jlambda, jcurve, /SILENT
READCOL,'/Users/covey/Dropbox/idl/data/twomass/twomass_H.par', FORMAT = '(X, F, F)', hlambda, hcurve, /SILENT
READCOL,'/Users/covey/Dropbox/idl/data/twomass/twomass_Ks.par', FORMAT = '(X, F, F)',klambda, kcurve, /SILENT

;read in WISE filter curves
READCOL,'/Users/covey/Dropbox/idl/data/WISE/RSR-W1.EE.txt', FORMAT = '(F, F)', W1lambda, W1curve, /SILENT
READCOL,'/Users/covey/Dropbox/idl/data/WISE/RSR-W2.EE.txt', FORMAT = '(F, F)', W2lambda, W2curve, /SILENT
READCOL,'/Users/covey/Dropbox/idl/data/WISE/RSR-W3.EE.txt', FORMAT = '(F, F)', W3lambda, W3curve, /SILENT
READCOL,'/Users/covey/Dropbox/idl/data/WISE/RSR-W4.EE.txt', FORMAT = '(F, F)', W4lambda, W4curve, /SILENT

;convert WISE wavelengths from microns to angstroms
W1lambda = W1lambda*10000.
W2lambda = W2lambda*10000.
W3lambda = W3lambda*10000.
W4lambda = W4lambda*10000.

END

;------------------------------------------------------------------------
FUNCTION jansk2vegamag, janskys

;calculate vega mags for UBVRI using fluxes in janskys (with
;zeropoints as given on
;http://www.astro.utoronto.ca/~patton/astro/mags.html
;reportedly using data from Bessel 1979)

Umag = -2.5*ALOG10(janskys[0]/1810.)
Bmag = -2.5*ALOG10(janskys[1]/4260.)
Vmag = -2.5*ALOG10(janskys[2]/3640.)
Rmag = -2.5*ALOG10(janskys[3]/3080.)
Imag = -2.5*ALOG10(janskys[4]/2550.)

mags=[Umag,Bmag,Vmag,Rmag,Imag]

RETURN, mags

END

;------------------------------------------------------------------------

function jansk2sdssmag, janskys
;sdss softening parameters
b_u = 1.4D-10
b_g = 0.9D-10
b_r = 1.2D-10
b_i = 1.8D-10
b_z = 7.4D-10

;calculate asinh AB mags

umag = -1.*(2.5/alog(10))*[asinh((janskys[0]/3631.)/(2*b_u))+alog(b_u)]
gmag = -1.*(2.5/alog(10))*[asinh((janskys[1]/3631.)/(2*b_g))+alog(b_g)]
rmag = -1.*(2.5/alog(10))*[asinh((janskys[2]/3631.)/(2*b_r))+alog(b_r)]
imag = -1.*(2.5/alog(10))*[asinh((janskys[3]/3631.)/(2*b_i))+alog(b_i)]
zmag = -1.*(2.5/alog(10))*[asinh((janskys[4]/3631.)/(2*b_z))+alog(b_z)]

mags=[umag,gmag,rmag,imag,zmag]

;SDSS is not quite on AB system.  Here are corrections from M. Blanton
;note that we are actually correcting our measurements from a
;'perfect' AB system onto the imperfect SDSS system since we used the 
;AB system 3631 zeropoint to get our flux densities originally.  So we
;apply our correction in the sense that m_AB-deltam = m_SDSS, opposite
;how you might naively think we should... (corrected by K. Covey
;6-13-05)

mags= mags-[-0.036, 0.012, 0.010, 0.028, 0.04]

return, mags

end

;----------------------------------------------------------------------

function response,ilambda, ivac_pass, lambda
;this interpolates over the SDSS filter curves and matches the number
;of elements in the data

response=SPLINE(ilambda, ivac_pass, lambda)

return,response

end

;----------------------------------------------------------------------

FUNCTION filter, lambda, flux, ulambda, ucurve, glambda, gcurve, rlambda, rcurve, ilambda, icurve, zlambda, zcurve, $
ubesslambda, ubesscurve, bbesslambda, bbesscurve, vbesslambda, vbesscurve, rbesslambda, rbesscurve, ibesslambda, ibesscurve, $
int_wfc_u_lambda, int_wfc_u_trans, int_wfc_g_lambda, int_wfc_g_trans, int_wfc_r_lambda, int_wfc_r_trans, $
int_wfc_i_lambda, int_wfc_i_trans, int_wfc_z_lambda, int_wfc_z_trans, $
Asahi_u_lambda, Asahi_u_trans, Asahi_g_lambda, Asahi_g_trans, Asahi_r_lambda, Asahi_r_trans, Asahi_i_lambda, Asahi_i_trans, Asahi_z_lambda, Asahi_z_trans, $
Andover_u_lambda, Andover_u_trans, Andover_g_lambda, Andover_g_trans, Andover_r_lambda, Andover_r_trans, Andover_i_lambda, Andover_i_trans, Andover_z_lambda, Andover_z_trans, Mould_r_lambda, Mould_r_trans

;<reading in filter info excised to separate function to save time;
;see specmags.pro for the original structure>

;these indices define where the data overlap the curves. For sdss
;spectra, only the g, r, and i filters completely.  If you are using
;spectra that span different wavelength ranges you will need to alter
;them.  I have commented out u and z indices in the ideal case that
;you span the entire wavelength range of sdss filters

;uind=WHERE(lambda LE MAX(ulambda))

;PRINT, Asahi_u_lambda
;PRINT, ulambda


uind = WHERE(lambda GE MIN(ulambda) AND lambda LE MAX(ulambda))
gind = WHERE(lambda GE MIN(glambda) AND lambda LE MAX(glambda))
rind = WHERE(lambda GE MIN(rlambda) AND lambda LE MAX(rlambda))
iind = WHERE(lambda GE MIN(ilambda) AND lambda LE max(ilambda))
zind = WHERE(lambda ge min(zlambda)and lambda le max(zlambda))

Asahi_uind = WHERE(lambda GE MIN(Asahi_u_lambda) AND lambda LE MAX(Asahi_u_lambda))
Asahi_gind = WHERE(lambda GE MIN(Asahi_g_lambda) AND lambda LE MAX(Asahi_g_lambda))
Asahi_rind = WHERE(lambda GE MIN(Asahi_r_lambda) AND lambda LE MAX(Asahi_r_lambda))
Asahi_iind = WHERE(lambda GE MIN(Asahi_i_lambda) AND lambda LE max(Asahi_i_lambda))
Asahi_zind = WHERE(lambda ge min(Asahi_z_lambda)and lambda le max(Asahi_z_lambda))

Andover_uind = WHERE(lambda GE MIN(Andover_u_lambda) AND lambda LE MAX(Andover_u_lambda))
Andover_gind = WHERE(lambda GE MIN(Andover_g_lambda) AND lambda LE MAX(Andover_g_lambda))
Andover_rind = WHERE(lambda GE MIN(Andover_r_lambda) AND lambda LE MAX(Andover_r_lambda))
Andover_iind = WHERE(lambda GE MIN(Andover_i_lambda) AND lambda LE max(Andover_i_lambda))
Andover_zind = WHERE(lambda ge min(Andover_z_lambda)and lambda le max(Andover_z_lambda))

int_wfc_uind = WHERE(lambda GE MIN(int_wfc_u_lambda) AND lambda LE MAX(int_wfc_u_lambda))
int_wfc_gind = WHERE(lambda GE MIN(int_wfc_g_lambda) AND lambda LE MAX(int_wfc_g_lambda))
int_wfc_rind = WHERE(lambda GE MIN(int_wfc_r_lambda) AND lambda LE MAX(int_wfc_r_lambda))
int_wfc_iind = WHERE(lambda GE MIN(int_wfc_i_lambda) AND lambda LE MAX(int_wfc_i_lambda))
int_wfc_zind = WHERE(lambda ge min(int_wfc_z_lambda) AND lambda LE MAX(int_wfc_z_lambda))

;zind = where(lambda ge min(zlambda))

Ubessind = WHERE(lambda GE MIN(Ubesslambda) AND lambda LE MAX(Ubesslambda))
Bbessind = WHERE(lambda GE MIN(Bbesslambda) AND lambda LE MAX(Bbesslambda))
Vbessind = WHERE(lambda GE MIN(Vbesslambda) AND lambda LE MAX(Vbesslambda))
Rbessind = WHERE(lambda GE MIN(Rbesslambda) AND lambda LE MAX(Rbesslambda))
Ibessind = WHERE(lambda GE MIN(Ibesslambda) AND lambda LE MAX(Ibesslambda))

Rmouldind = WHERE(lambda GE MIN(mould_r_lambda) AND lambda LE MAX(mould_r_lambda))


;;these are testing where the filters overlap the data (comment out in
;;the ideal case of having full wavelength coverage)
;uresind = WHERE(ulambda GE MIN(lambda) AND ulambda LE MAX(lambda))
;zresind = WHERE(zlambda LE MAX(lambda) AND zlambda GE MIN(lambda))
;gresind = WHERE(glambda GE MIN(lambda) AND glambda LE MAX(lambda))
;Ubessresind = WHERE(Ubesslambda GE MIN(lambda))

;calculate the response for each band
ures = RESPONSE(ulambda, ucurve, lambda[uind])
gres = RESPONSE(glambda, gcurve, lambda[gind])
rres = RESPONSE(rlambda, rcurve, lambda[rind])
ires = RESPONSE(ilambda, icurve, lambda[iind])
zres = RESPONSE(zlambda, zcurve, lambda[zind])

Asahi_ures = RESPONSE(Asahi_u_lambda, Asahi_u_trans, lambda[Asahi_uind])
Asahi_gres = RESPONSE(Asahi_g_lambda, Asahi_g_trans, lambda[Asahi_gind])
Asahi_rres = RESPONSE(Asahi_r_lambda, Asahi_r_trans, lambda[Asahi_rind])
Asahi_ires = RESPONSE(Asahi_i_lambda, Asahi_i_trans, lambda[Asahi_iind])
Asahi_zres = RESPONSE(Asahi_z_lambda, Asahi_z_trans, lambda[Asahi_zind])

Andover_ures = RESPONSE(Andover_u_lambda, Andover_u_trans, lambda[Andover_uind])
Andover_gres = RESPONSE(Andover_g_lambda, Andover_g_trans, lambda[Andover_gind])
Andover_rres = RESPONSE(Andover_r_lambda, Andover_r_trans, lambda[Andover_rind])
Andover_ires = RESPONSE(Andover_i_lambda, Andover_i_trans, lambda[Andover_iind])
Andover_zres = RESPONSE(Andover_z_lambda, Andover_z_trans, lambda[Andover_zind])

int_wfc_ures = RESPONSE(int_wfc_u_lambda, int_wfc_u_trans, lambda[int_wfc_uind])
int_wfc_gres = RESPONSE(int_wfc_g_lambda, int_wfc_g_trans, lambda[int_wfc_gind])
int_wfc_rres = RESPONSE(int_wfc_r_lambda, int_wfc_r_trans, lambda[int_wfc_rind])
int_wfc_ires = RESPONSE(int_wfc_i_lambda, int_wfc_i_trans, lambda[int_wfc_iind])
int_wfc_zres = RESPONSE(int_wfc_z_lambda, int_wfc_z_trans, lambda[int_wfc_zind])


Ubessres = RESPONSE(Ubesslambda, Ubesscurve, lambda[Ubessind])
Bbessres = RESPONSE(Bbesslambda, Bbesscurve, lambda[Bbessind])
Vbessres = RESPONSE(Vbesslambda, Vbesscurve, lambda[Vbessind])
Rbessres = RESPONSE(Rbesslambda, Rbesscurve, lambda[Rbessind])
Ibessres = RESPONSE(Ibesslambda, Ibesscurve, lambda[Ibessind])

Rmouldres = RESPONSE(mould_r_lambda, mould_r_trans, lambda[rmouldind])


;take a flux at each wavelength in units of ergs/cm^2-sec-Ang,
;convolve with filter response, use the fact that dn=dl^2/c where dl
;is a change in wavelength and dn is a change in frequency to convert
;the 1/A to 1/Hz (remember that 2.998d18 is c in angstroms per second
;to make the units work out right), then use the fact that 1 Jy =
;10^-23 ergs/sec/cm^2/Hz to get the flux at each wavelength in janskys.

uflux = DOUBLE(ures*flux[uind]*(lambda[uind]^2/2.998d18)*1.d23)
gflux = DOUBLE(gres*flux[gind]*(lambda[gind]^2/2.998d18)*1.d23)
rflux = DOUBLE(rres*flux[rind]*(lambda[rind]^2/2.998d18)*1.d23)
iflux = DOUBLE(ires*flux[iind]*(lambda[iind]^2/2.998d18)*1.d23)
zflux = DOUBLE(zres*flux[zind]*(lambda[zind]^2/2.998d18)*1.d23)

Asahi_uflux = DOUBLE(Asahi_ures*flux[Asahi_uind]*(lambda[Asahi_uind]^2/2.998d18)*1.d23)
Asahi_gflux = DOUBLE(Asahi_gres*flux[Asahi_gind]*(lambda[Asahi_gind]^2/2.998d18)*1.d23)
Asahi_rflux = DOUBLE(Asahi_rres*flux[Asahi_rind]*(lambda[Asahi_rind]^2/2.998d18)*1.d23)
Asahi_iflux = DOUBLE(Asahi_ires*flux[Asahi_iind]*(lambda[Asahi_iind]^2/2.998d18)*1.d23)
Asahi_zflux = DOUBLE(Asahi_zres*flux[Asahi_zind]*(lambda[Asahi_zind]^2/2.998d18)*1.d23)

Andover_uflux = DOUBLE(Andover_ures*flux[Andover_uind]*(lambda[Andover_uind]^2/2.998d18)*1.d23)
Andover_gflux = DOUBLE(Andover_gres*flux[Andover_gind]*(lambda[Andover_gind]^2/2.998d18)*1.d23)
Andover_rflux = DOUBLE(Andover_rres*flux[Andover_rind]*(lambda[Andover_rind]^2/2.998d18)*1.d23)
Andover_iflux = DOUBLE(Andover_ires*flux[Andover_iind]*(lambda[Andover_iind]^2/2.998d18)*1.d23)
Andover_zflux = DOUBLE(Andover_zres*flux[Andover_zind]*(lambda[Andover_zind]^2/2.998d18)*1.d23)

int_wfc_uflux = DOUBLE(int_wfc_ures*flux[int_wfc_uind]*(lambda[int_wfc_uind]^2/2.998d18)*1.d23)
int_wfc_gflux = DOUBLE(int_wfc_gres*flux[int_wfc_gind]*(lambda[int_wfc_gind]^2/2.998d18)*1.d23)
int_wfc_rflux = DOUBLE(int_wfc_rres*flux[int_wfc_rind]*(lambda[int_wfc_rind]^2/2.998d18)*1.d23)
int_wfc_iflux = DOUBLE(int_wfc_ires*flux[int_wfc_iind]*(lambda[int_wfc_iind]^2/2.998d18)*1.d23)
int_wfc_zflux = DOUBLE(int_wfc_zres*flux[int_wfc_zind]*(lambda[int_wfc_zind]^2/2.998d18)*1.d23)


Ubessflux = DOUBLE(Ubessres*flux[Ubessind]*(lambda[Ubessind]^2/2.998d18)*1.d23)
Bbessflux = DOUBLE(Bbessres*flux[Bbessind]*(lambda[Bbessind]^2/2.998d18)*1.d23)
Vbessflux = DOUBLE(Vbessres*flux[Vbessind]*(lambda[Vbessind]^2/2.998d18)*1.d23)
Rbessflux = DOUBLE(Rbessres*flux[Rbessind]*(lambda[Rbessind]^2/2.998d18)*1.d23)
Ibessflux = DOUBLE(Ibessres*flux[Ibessind]*(lambda[Ibessind]^2/2.998d18)*1.d23)

Rmouldflux = DOUBLE(Rmouldres*flux[Rmouldind]*(lambda[Rmouldind]^2/2.998d18)*1.d23)


;integrate using trapezoidal rule to get the total flux (change 
;the ordinate into frequency space by using c [in angstroms per 
;second] / wavelength to get Hz) and then divide by the "effective" 
;filter width (which is essentially the width a perfect [response = 1]
;boxcar filter would have to have to get the same amount of flux
;through as our real filter did) to get the flux density in Janskys 

umag = ABS(TSUM(2.998d18/lambda[uind],uflux)/(2.998d18*TSUM(lambda[uind],ures)/3580.^2))
gmag = ABS(TSUM(2.998d18/lambda[gind],gflux)/(2.998d18*TSUM(lambda[gind],gres)/4755.^2))
rmag = ABS(TSUM(2.998d18/lambda[rind],rflux)/(2.998d18*TSUM(lambda[rind],rres)/6205.^2)) ;6205
imag = ABS(TSUM(2.998d18/lambda[iind],iflux)/(2.998d18*TSUM(lambda[iind],ires)/7480.^2))
zmag = ABS(TSUM(2.998d18/lambda[zind],zflux)/(2.998d18*TSUM(lambda[zind],zres)/8855.^2))

Asahi_umag = ABS(TSUM(2.998d18/lambda[Asahi_uind],Asahi_uflux)/(2.998d18*TSUM(lambda[Asahi_uind],Asahi_ures)/3580.^2))
Asahi_gmag = ABS(TSUM(2.998d18/lambda[Asahi_gind],Asahi_gflux)/(2.998d18*TSUM(lambda[Asahi_gind],Asahi_gres)/4755.^2))
Asahi_rmag = ABS(TSUM(2.998d18/lambda[Asahi_rind],Asahi_rflux)/(2.998d18*TSUM(lambda[Asahi_rind],Asahi_rres)/6205.^2)) ;6205
Asahi_imag = ABS(TSUM(2.998d18/lambda[Asahi_iind],Asahi_iflux)/(2.998d18*TSUM(lambda[Asahi_iind],Asahi_ires)/7480.^2))
Asahi_zmag = ABS(TSUM(2.998d18/lambda[Asahi_zind],Asahi_zflux)/(2.998d18*TSUM(lambda[Asahi_zind],Asahi_zres)/8855.^2))

Andover_umag = ABS(TSUM(2.998d18/lambda[Andover_uind],Andover_uflux)/(2.998d18*TSUM(lambda[Andover_uind],Andover_ures)/3580.^2))
Andover_gmag = ABS(TSUM(2.998d18/lambda[Andover_gind],Andover_gflux)/(2.998d18*TSUM(lambda[Andover_gind],Andover_gres)/4755.^2))
Andover_rmag = ABS(TSUM(2.998d18/lambda[Andover_rind],Andover_rflux)/(2.998d18*TSUM(lambda[Andover_rind],Andover_rres)/6205.^2)) ;6205
Andover_imag = ABS(TSUM(2.998d18/lambda[Andover_iind],Andover_iflux)/(2.998d18*TSUM(lambda[Andover_iind],Andover_ires)/7480.^2))
Andover_zmag = ABS(TSUM(2.998d18/lambda[Andover_zind],Andover_zflux)/(2.998d18*TSUM(lambda[Andover_zind],Andover_zres)/8855.^2))

int_wfc_umag = ABS(TSUM(2.998d18/lambda[int_wfc_uind],int_wfc_uflux)/(2.998d18*TSUM(lambda[int_wfc_uind],int_wfc_ures)/3580.^2))
int_wfc_gmag = ABS(TSUM(2.998d18/lambda[int_wfc_gind],int_wfc_gflux)/(2.998d18*TSUM(lambda[int_wfc_gind],int_wfc_gres)/4755.^2))
int_wfc_rmag = ABS(TSUM(2.998d18/lambda[int_wfc_rind],int_wfc_rflux)/(2.998d18*TSUM(lambda[int_wfc_rind],int_wfc_rres)/6205.^2)) ;6205
int_wfc_imag = ABS(TSUM(2.998d18/lambda[int_wfc_iind],int_wfc_iflux)/(2.998d18*TSUM(lambda[int_wfc_iind],int_wfc_ires)/7480.^2))
int_wfc_zmag = ABS(TSUM(2.998d18/lambda[int_wfc_zind],int_wfc_zflux)/(2.998d18*TSUM(lambda[int_wfc_zind],int_wfc_zres)/8855.^2))

;try this calculation to see if it works better.  NOPE.  It DOESNT.
;umag = ABS(TSUM(3.d18/lambda[uind],uflux)/TSUM(3.d18/lambda[uind],ures))
;gmag = ABS(TSUM(3.d18/lambda[gind],gflux)/TSUM(3.d18/lambda[gind],gres))
;rmag = ABS(TSUM(3.d18/lambda[rind],rflux)/TSUM(3.d18/lambda[rind],rres)) ;6205
;imag = ABS(TSUM(3.d18/lambda[iind],iflux)/TSUM(3.d18/lambda[iind],ires))
;zmag = ABS(TSUM(3.d18/lambda[zind],zflux)/TSUM(3.d18/lambda[zind],zres))

;SET_PLOT, 'X'

;PRINT, (2.998d18*TSUM(lambda[iind],ires)/7480.^2), TSUM(3.d18/lambda[iind],ires)
;wait = get_kbrd(1)
;PLOT, lambda[iind],ires
;wait = get_kbrd(1)
;PLOT, 3.d18/lambda[iind],ires

Ubessmag = ABS(TSUM(3.d18/lambda[Ubessind],Ubessflux)/(2.998d18*TSUM(lambda[Ubessind],Ubessres)/3600.^2))
Bbessmag = ABS(TSUM(3.d18/lambda[Bbessind],Bbessflux)/(2.998d18*TSUM(lambda[Bbessind],Bbessres)/4400.^2))
Vbessmag = ABS(TSUM(3.d18/lambda[Vbessind],Vbessflux)/(2.998d18*TSUM(lambda[Vbessind],Vbessres)/5500.^2))
Rbessmag = ABS(TSUM(3.d18/lambda[Rbessind],Rbessflux)/(2.998d18*TSUM(lambda[Rbessind],Rbessres)/6500.^2))
Ibessmag = ABS(TSUM(3.d18/lambda[Ibessind],Ibessflux)/(2.998d18*TSUM(lambda[Ibessind],Ibessres)/8100.^2))

Rmouldmag = ABS(TSUM(3.d18/lambda[Rmouldind],Rmouldflux)/(2.998d18*TSUM(lambda[Rmouldind],Rmouldres)/6500.^2))


fluxes = FLTARR(6,5)
fluxes[0,*] = [umag,gmag,rmag,imag,zmag]
fluxes[1,*] = [Ubessmag,Bbessmag,Vbessmag,Rbessmag,Ibessmag]
fluxes[2,*] = [Asahi_umag,Asahi_gmag,Asahi_rmag,Asahi_imag,Asahi_zmag]
fluxes[3,*] = [Andover_umag,Andover_gmag,Andover_rmag,Andover_imag,Andover_zmag]
fluxes[4,*] = [int_wfc_umag,int_wfc_gmag,int_wfc_rmag,int_wfc_imag,int_wfc_zmag]
fluxes[5,3] = Rmouldmag

;PRINT, fluxes

;fluxes 
RETURN,fluxes

END

;---------------------------------------------------------------------

FUNCTION READ_PHOENIX_SPECTRUM, filename

spectrum = READFITS(filename, header, /SILENT)

n_spectrum = N_ELEMENTS(spectrum)

pixel_zero = SXPAR(header, 'CRVAL1')
pixel_step = SXPAR(header, 'CDELT1')

for_lambda = FINDGEN(n_spectrum)

lambda = EXP(pixel_zero+pixel_step*for_lambda)

datastruct = DBLARR(n_spectrum,2)
datastruct[*,0] = lambda
datastruct[*,1] = spectrum

RETURN, datastruct

END

;;-------

PRO downsample_Phoenix_specs, new_res

PRINT, 'USING THE RIGHT CODE'
wait = GET_KBRD(1)

LOAD_ALL_FILTERS, ulambda, ucurve, glambda, gcurve, rlambda, rcurve, ilambda, icurve, zlambda, zcurve, $
ubesslambda, ubesscurve, bbesslambda, bbesscurve, vbesslambda, vbesscurve, rbesslambda, rbesscurve, ibesslambda, ibesscurve, $
int_wfc_u_lambda, int_wfc_u_trans, int_wfc_g_lambda, int_wfc_g_trans, int_wfc_r_lambda, int_wfc_r_trans, $
int_wfc_i_lambda, int_wfc_i_trans, int_wfc_z_lambda, int_wfc_z_trans, $
Asahi_u_lambda, Asahi_u_trans, Asahi_g_lambda, Asahi_g_trans, Asahi_r_lambda, Asahi_r_trans, Asahi_i_lambda, Asahi_i_trans, Asahi_z_lambda, Asahi_z_trans, $
Andover_u_lambda, Andover_u_trans, Andover_g_lambda, Andover_g_trans, Andover_r_lambda, Andover_r_trans, Andover_i_lambda, Andover_i_trans, Andover_z_lambda, Andover_z_trans, jlambda, jcurve, hlambda, hcurve, klambda, kcurve, W1lambda, W1curve, W2lambda, W2curve, W3lambda, W3curve, W4lambda, W4curve, Mould_r_lambda, Mould_r_trans

fake_spec = DBLARR(37132) ;-- for R = 3500.

;fake_spec = DBLARR(317608) ; -- for R = 30000.
fake_spec[*] = -99.

structdef = {PhoenixName: '', teff: -99., feh: -99., logg: -99., $
             lambda: fake_spec, flux: fake_spec, resolution: new_res, $
             sdss_u: -99., sdss_g: -99., sdss_r: -99., sdss_i: -99., sdss_z: -99., $
             int_wfc_u: -99., int_wfc_g: -99., int_wfc_r: -99., int_wfc_i: -99.,int_wfc_z: -99., $ 
             two_j: -99., two_h: -99., two_k: -99., W1: -99., W2: -99., W3: -99., W4: -99., $
             Bessell_U: -99., Bessell_B: -99., Bessell_V: -99., Bessell_R: -99., Bessell_I: -99.}
;model_struct = REPLICATE(structdef, 5343)
model_struct = REPLICATE(structdef, 1)

temps = (FINDGEN(47)+23)*100
;temps = (FINDGEN(47)+35)*100
hotter = (FINDGEN(26)*200)+7000

all_temps = [temps,hotter]

n_temps = N_ELEMENTS(all_temps)

feh = [-2., -1.5, -1.0, -0.5, -0, 0.5, 1.0]

models = 0

FOR i=0,n_temps-1 DO BEGIN

   teff_string = STRTRIM(STRING(ROUND(all_temps[i])),2)
   ;PRINT, STRLEN(teff_string)
   IF STRLEN(teff_string) EQ 4 THEN teff_string = '0'+teff_string

   PRINT, teff_string

   FOR j=0, 6 DO BEGIN

      IF feh[j] LT 0 THEN feh_string = STRMID(STRTRIM(STRING(feh[j]),2),0,4) ELSE IF feh[j] GT 0 THEN feh_string = '+'+STRMID(STRTRIM(STRING(feh[j]),2),0,3) ELSE IF feh[j] EQ 0 THEN feh_string = '-'+STRMID(STRTRIM(STRING(feh[j]),2),0,3)

      FOR k=0,11 DO BEGIN

         logg = k*0.5

        filename = 'lte'+teff_string+'-'+STRMID(STRTRIM(STRING(logg),2),0,4)+feh_string+'.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'

        ;PRINT, filename

        ;should test if the model exists first!
        IF FILE_TEST(filename) EQ 1 THEN BEGIN
           ;PRINT, 'found it!' 
           model_struct.PhoenixName = filename
           model_struct.teff = all_temps[i]
           model_struct.feh = feh[j]
           model_struct.logg = k*0.5

           ;PRINT, model_struct[models].logg

           this_model = READ_PHOENIX_SPECTRUM(filename)

           downsampled_model = downsample(this_model[*,0], this_model[*,1], new_res) 
           model_struct.lambda = downsampled_model[*,0]
           model_struct.flux = downsampled_model[*,1]
           PRINT, N_ELEMENTS(this_model[*,0])
           PRINT, N_ELEMENTS(downsampled_model[*,0])

           optical_phot = FILTER(this_model[*,0], this_model[*,1]*1E-8, ulambda, ucurve, glambda, gcurve, rlambda, rcurve, ilambda, icurve, zlambda, zcurve, $
ubesslambda, ubesscurve, bbesslambda, bbesscurve, vbesslambda, vbesscurve, rbesslambda, rbesscurve, ibesslambda, ibesscurve, $
int_wfc_u_lambda, int_wfc_u_trans, int_wfc_g_lambda, int_wfc_g_trans, int_wfc_r_lambda, int_wfc_r_trans, $
int_wfc_i_lambda, int_wfc_i_trans, int_wfc_z_lambda, int_wfc_z_trans, $
Asahi_u_lambda, Asahi_u_trans, Asahi_g_lambda, Asahi_g_trans, Asahi_r_lambda, Asahi_r_trans, Asahi_i_lambda, Asahi_i_trans, Asahi_z_lambda, Asahi_z_trans, $
Andover_u_lambda, Andover_u_trans, Andover_g_lambda, Andover_g_trans, Andover_r_lambda, Andover_r_trans, Andover_i_lambda, Andover_i_trans, Andover_z_lambda, Andover_z_trans, Mould_r_lambda, Mould_r_trans)
           sdss_mags = jansk2sdssmag(optical_phot[0,*])
           model_struct.sdss_u = sdss_mags[0]
           model_struct.sdss_g = sdss_mags[1]
           model_struct.sdss_r = sdss_mags[2]
           model_struct.sdss_i = sdss_mags[3]
           model_struct.sdss_z = sdss_mags[4]

           int_wfc_mags = jansk2sdssmag(optical_phot[4,*])
           model_struct.int_wfc_u = int_wfc_mags[0]
           model_struct.int_wfc_g = int_wfc_mags[1]
           model_struct.int_wfc_r = int_wfc_mags[2]
           model_struct.int_wfc_i = int_wfc_mags[3]
           model_struct.int_wfc_z = int_wfc_mags[4]

           IR_phot = IRFILTER(this_model[*,0], this_model[*,1]*1E-8, jlambda, jcurve, hlambda, hcurve, klambda, kcurve, W1lambda, W1curve, W2lambda, W2curve, W3lambda, W3curve, W4lambda, W4curve)
           IR_mags = jansk2IRmag(IR_phot)
           model_struct.two_j = IR_mags[0]
           model_struct.two_h = IR_mags[1]
           model_struct.two_k = IR_mags[2]
;;           model_struct[i].W1 = IR_mags[3]
;;           model_struct[i].W2 = IR_mags[4]
;;           model_struct[i].W3 = IR_mags[5]
;;           model_struct[i].W4 = IR_mags[6]

           ;PRINT, model_struct[i].sdss_i -  model_struct[i].sdss_z, model_struct[i].sdss_r -  model_struct[i].sdss_i

           UBVRI_mags = jansk2vegamag(optical_phot[1,*])
           model_struct.Bessell_U = UBVRI_mags[0]
           model_struct.Bessell_B = UBVRI_mags[1]
           model_struct.Bessell_V = UBVRI_mags[2]
           model_struct.Bessell_R = UBVRI_mags[3]
           model_struct.Bessell_I = UBVRI_mags[4]
 
           MWRFITS, model_struct, 'lte'+teff_string+'-'+STRMID(STRTRIM(STRING(logg),2),0,4)+feh_string+'.PHOENIX-ACES-AGSS-COND-2011-R-'+STRMID(STRTRIM(STRING(new_res),2),0,5)+'fits', /CREATE

;           models = models+1

        ENDIF 

       ENDFOR

   ENDFOR

ENDFOR

END
