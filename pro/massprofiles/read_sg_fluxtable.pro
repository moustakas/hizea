function read_sg_fluxtable
; jm17jul03siena - read the table produced by Aleks' student, Sophia, and
; pack it into a format that massprofile_isedfit can use.

    massdir = massprofiles_path()
    parent = mrdfits(massdir+'massprofiles_photometry.fits.gz',1)
    
    filterlist = massprofiles_filterlist()
    nfilt = n_elements(filterlist)
    
    fluxfile = getenv('HIZEA_DIR')+'/etc/sg_fluxtable_nm.txt'
    readcol, fluxfile, id, f475, f475_ivar, f814, f814_ivar, f1600, f1600_ivar, z, $
      format='A,F,F,F,F,F,F,D', /silent
    nflux = n_elements(f475)

    gal = strmid(id,0,5)
    galaxy = gal[uniq(gal,sort(gal))]
    ap = strmid(id,6,2)

    cat = {$
      galaxy: '', $
      ra: 0.0D, $
      dec: 0.0D, $
      z: 0.0D, $
      aperture: '', $
      maggies: fltarr(nfilt), $
      ivarmaggies: fltarr(nfilt)}
    cat = replicate(cat,nflux)

;   cat.galaxy = gal
;   cat.z = z
    cat.aperture = ap
    cat.maggies = 1E-9*transpose([[f475*0.0], [f475], [f814], [f1600]])
    cat.ivarmaggies = 1E18*transpose([[f475_ivar*0.0], [f475_ivar], [f814_ivar], [f1600_ivar]])

; match back to the parent sample
    for ii = 0, n_elements(parent)-1 do begin
       ww = where(strmid(parent[ii].galaxy,0,5) eq gal)
       cat[ww].galaxy = parent[ii].galaxy
       cat[ww].ra = parent[ii].ra
       cat[ww].dec = parent[ii].dec
       cat[ww].z = parent[ii].z
    endfor

; sort by redshift (for iSEDfit)
    cat = cat[sort(cat.z)]
    
return, cat
end
    
