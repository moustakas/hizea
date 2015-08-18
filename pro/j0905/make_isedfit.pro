FUNCTION get_data, field
;  ;; Start with the 0d catalog, and grab the minimal set of data.
;  ;; You will want to adjust this to include the right 
;  data = primus_read_zerod(field=field, rerun=2165, phot=photo)
;  
;  ii = where(data.zprimus GT 0.2 AND $
;             data.zprimus LE 1.2 AND $
;             data.zprimus_zconf GE 3, nii)
;  
;  ;; Here is a pretty quick hack that just removes spitzer.  You should
;  ;; check to see if these are the filters that you want to use
;  ;; and implement a better technique to get the right filters for
;  ;; each field. Right now it probably is not.
;  filters = photo[0].filterlist
;  jj = where(strmatch(filters, '*spitzer*') EQ 0, njj)
;  
;  splog, 'Using data from these filters: ', strjoin(filters[jj], ', ')
;  
;  
;  ;; Create a new structure for the ones that we want
;  cat = {ra:-1D, $
;         dec:-1D, $
;         z: -1.0, $
;         ;; you probably want to add more tags here so that you can
;         ;; use this data structure without always matching and dealing
;         ;; with the slow primus_read_zerod function.
;         filters:filters, $
;         maggies:fltarr(njj)-1, $
;         maggiesivar:fltarr(njj)-1 }
;  cat = replicate(cat, nii)
;  cat.ra = data[ii].ra
;  cat.dec = data[ii].dec
;  cat.z = data[ii].zprimus
;  cat.maggies = photo[ii].maggies[jj]
;  cat.ivarmaggies = photo[ii].ivarmaggies[jj]
;  return,  cat

; amd130717 -- adding in data for foreground galaxy
; note mendez is using r807 of impro
ra = 136.35035122d
dec = 57.98775331d

mag = [21.979784, 21.743961, 20.45031, 20.062166, 19.538656]
mag_err = [0.322197, 0.106713, 0.052093, 0.053842, 0.163799]
redshift = 0.4134

maggies = 10.^(mag*(-1.)/2.5)
maggies_err = mag_err * maggies / 1.086
maggies_ivar = 1. / (maggies_err)^2

filterlist = sdss_filterlist()

  cat = {ra: ra, $
         dec: dec, $
         z: redshift, $
         filters: filterlist, $
         maggies: maggies, $
         maggiesivar: maggies_ivar }

return, cat
END





PRO make_isedfit, prelim=prelim, $
                  models=models, $
                  montegrid=montegrid, $
                  isedfit=isedfit, $
                  qaplot=qaplot, $
                  clobber=clobber
  ;prefix = 'cosmos'
  ;ver = 'v5.2'
  prefix = 'fg'  

  ;isedfit_dir = '/path/to/an/existing/directory'
  isedfit_dir = '/Users/aleks/hizea/p2/isedfit/'
  isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
  sfhgrid_paramfile = isedfit_dir+prefix+'_sfhgrid.par'
  supergrid_paramfile = isedfit_dir+prefix+'_supergrid.par'
  
; get the data
  data = get_data(prefix)
  cat = data  

; redshift bins
  zmin = 0.01
  zmax = 2.5
  nzz = 100
  
; SFHgrid priors                  [Previous PRIMUS settings]
  nage = 20                     ; 50
  nmonte = 1000                 ; 1000
  Z = [0.004, 0.04]             ; [0.004, 0.03]
  minage = 0.1                  ; 0.1
  maxage = 13.0                 ; 13.0
  flatAV = 0                    ; 0
  AV = [0.35, 2.0]              ; [0.35 2.0]
  tau = [0.01, 5.0]             ; [0.01,1.0]
  delayed = 1                   ; 1
  preset_bursts = 1             ; 1

; amd130717 -- mendez says "you might want to use delayed=0,
;              although only affects at z>2"
  
; supergrid parameters    
  sfhgrid = 1
  supergrid = 1
  synthmodels = 'fsps'
  
  imf = 'chab'
  redcurve = 1                  ; Charlot & Fall
  
; output filename reconstruction
  outfile = isedfit_dir+$
            prefix+'_'+$
            synthmodels+'_'+$
            imf+'_charlot_sfhgrid'+$
            strtrim(string(supergrid, format='(I02)'), 2)+$
            '.fits.gz'
  
; --------------------------------------------------
; do the preliminaries: build the parameter files:: 
; iSEDfit requires three different parameter files. A global parameter
; file that is tied to a particular dataset (e.g., it has filter
; information); a ”sfhgrid” (star formation history grid) parameter
; file, which specifies the prior parameter choices on the star
; formation history, metallicity, attenuation, burst parameters,
; etc.; and a ”supergrid” parameter file, which specifies the SPS
; models, IMF, and attenuation/reddening curve.

  IF keyword_set(prelim) THEN BEGIN
    ;; This makes the global parameter file. Choose the filter list, 
    ;; and be sure they have been incorporated into K-correct. 
    ;; Choose the minimum and maximum redshifts you will consider. 
    ;; Choose a unique prefix for the project.
    write_isedfit_paramfile, cat[0].filters, prefix=prefix, minz=zmin, $
       maxz=zmax, nzz=nzz, igm=1, isedfit_dir=isedfit_dir, clobber=clobber
    
    
    ;; Next build the ”sfhgrid” parameter file. This file specifies your
    ;; choice of star formation history (SFH) priors. You can have as
    ;; many combinations of priors as you want and they will each be
    ;; assigned a unique identifier. 
    write_sfhgrid_paramfile, sfhgrid_paramfile, params, clobber=clobber, $
       sfhgrid=sfhgrid, nage=nage, nmonte=nmonte, Z=Z, minage=minage, $
       preset_bursts=preset_bursts, $
       maxage=maxage, AV=AV, tau=tau, flatAV=flatAV, delayed=delayed
    
    ;; Next build the ”supergrid” parameter file. This file relates each
    ;; SFHGRID to a particular choice of SPS models, IMF, and attenuation
    ;; curve. The strength of this data model is that you can easily
    ;; explore the effect of different SPS models or attenuation curve
    ;; on your results.
    write_supergrid_paramfile, supergrid_paramfile, supergrid=supergrid, $
       sfhgrid=sfhgrid, synthmodels=synthmodels, imf=imf, redcurve=redcurve, $
       clobber=clobber
  ENDIF
  
; --------------------------------------------------
; Build the monte carlo grid.
; The Monte Carlo based star formation history grids are used to sample
; the possible parameter spaces depending on the sfhgrid parameter file
; and SSPs.
;  [This takes some time, so start it and go to tea time.]
; amd130717 -- this took 20 minutes for sdss_filterlist()
  IF keyword_set(montegrid) THEN BEGIN
    build_montegrids,  sfhgrid_paramfile,  $
      supergrid_paramfile=supergrid_paramfile,  $
      isedfit_dir=isedfit_dir, $
      clobber=clobber
  ENDIF
    
  
  
  
; --------------------------------------------------
; build the models
; Synthesize photometry on a grid of redshift for all the
; models generated by build_montegrids.  Uses the iSEDfit
; parameter file which specifies the filter, the redshift
; range, the cosmological model (which affects the
; luminosity distance), and which specifies whether or
; not to include IGM attenuation.
; amd130717 -- this took 13 minutes
  IF keyword_set(models) THEN BEGIN
    isedfit_models, isedfit_paramfile, isedfit_dir=isedfit_dir, $
                    supergrid_paramfile=supergrid_paramfile, clobber=clobber
  ENDIF

; --------------------------------------------------
; do the fitting
; Model the spectral energy distributions of galaxies using
; population synthesis models to infer their physical properties. 
  IF keyword_set(isedfit) THEN BEGIN
    isedfit, isedfit_paramfile, cat.maggies, cat.maggiesivar, $;cat.ivarmaggies, $
      cat.z, result, isedfit_dir=isedfit_dir, clobber=clobber, $
      supergrid_paramfile=supergrid_paramfile, outprefix=outprefix, $
      sfhgrid_paramfile=sfhgrid_paramfile, galchunksize=1000L
  ENDIF

; --------------------------------------------------
; make a QAplot
; Generate quality-assurance plots from ISEDFIT output.
  IF keyword_set(qaplot) THEN  BEGIN
    ;ii = lindgen(n_elements(cat)) ;; possibly plot all of them
    ;index = ii[shuffle_indx(nii,num=50)];; randomly select 50 of them
    ;galaxy = 'index:'+number_formatter(index)
    
    index = 0L
    galaxy = 'fg'
 
    isedfit_qaplot, isedfit_paramfile, $
      supergrid_paramfile=supergrid_paramfile, $
      isedfit_dir=isedfit_dir, $
      montegrids_dir=montegrids_dir, $
      index=index, $
      galaxy=galaxy, $
      outprefix=outprefix, clobber=clobber, $
      /xlog
  ENDIF
    
END
