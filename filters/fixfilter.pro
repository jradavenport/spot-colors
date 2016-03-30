pro fixfilter

    ; need to read in all the filters, then spit out w/ same formatting

    sdss = ['u.dat.txt','g.dat.txt','r.dat.txt','i.dat.txt','z.dat.txt']
    kep = 'kepler_response_hires1'
    nir = ['sec6_4a.tbl1','sec6_4a.tbl2','sec6_4a.tbl3']

    out = ['u','g','r','i','z']
    for k=0,n_elements(sdss)-1 do begin
        ; SDSS filters are in the right wavelength, no resample needed
        readcol, sdss[k], wave, res, f='(A,X,X,A)', comm='#',/silent
        forprint, textout=out[k]+'_filter.txt', wave,res, /silent,/ nocomment, format='(F10.4, F10.5)'
    endfor



        out = ['j','h','k']
        for k=0,n_elements(nir)-1 do begin

            readcol, nir[k], wave, res, f='(F,F)',/silent
            forprint, textout=out[k]+'_filter.txt', wave * 1d4,res, /silent,/ nocomment, format='(F10.4, F10.5)'
        endfor

        readcol, kep, wave, res, comm='#',/silent
        forprint, textout='kep_filter.txt', wave * 1d1, res, /silent,/ nocomment, format='(F10.4, F10.5)'

    stop
end
