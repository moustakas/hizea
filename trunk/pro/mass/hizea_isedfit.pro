pro hizea_isedfit, models=models, isedfit=isedfit, qaplot=qaplot, $
  measure=measure, clobber=clobber, sdss=sdss, debug=debug, noirac=noirac
; jm10dec20ucsd - derive stellar masses for the HIZEA sample

    iopath = getenv('HIZEA_DATA')+'/analysis/mass/'
    paramfile = iopath+'hizea_isedfit.par'
    params = read_isedfit_paramfile(paramfile)
    
; --------------------------------------------------
; build the models
    if keyword_set(models) then isedfit_models, $
      paramfile, iopath=iopath, clobber=clobber

; --------------------------------------------------
; do the fitting!  
    if keyword_set(isedfit) then begin
       maggies = read_hizea_photometry(ivarmaggies=ivarmaggies,$
         zobj=zobj,filterlist=filt)
       if keyword_set(noirac) then begin
          toss = where(strmatch(filt,'*irac*'))
          ivarmaggies[toss,*] = 0.0
          outprefix = 'hizea_noirac'
       endif
       isedfit, paramfile, maggies, ivarmaggies, zobj, result, $
         iopath=iopath, outprefix=outprefix, galchunksize=galchunk, $
         clobber=clobber, debug=debug, index=index
    endif 

; --------------------------------------------------
; make a QAplot
    if keyword_set(qaplot) then begin
       if keyword_set(noirac) then outprefix = 'hizea_noirac'
       junk = read_hizea_photometry(galaxy=galaxy)
       isedfit_qaplot, paramfile, isedfit, iopath=iopath, galaxy=galaxy, $
         index=index, clobber=clobber, outprefix=outprefix
    endif

; --------------------------------------------------
; measure rest-frame quantities
    if keyword_set(measure) then begin
       isedfit_measure, paramfile, measure, isedfit, iopath=iopath, $
         clobber=clobber, outprefix=outprefix
    endif

return
end
