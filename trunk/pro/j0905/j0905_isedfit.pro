pro j0905_write_paramfile, paramfile, prefix=prefix, super=super

    zminmax = [0.71,0.72]
    nzz = 3
    igm = 1
    zlog = 0

    h100 = 0.7
    omega0 = 0.3
    omegal = 0.7

    filters = j0905_filterlist()
    sfhgridstring = strtrim(super.sfhgrid,2)
    redcurvestring = redcurve2string(super.redcurve)
    synthmodels = strtrim(super.synthmodels,2)
    imf = strtrim(super.imf,2)

    splog, 'Writing '+paramfile
    zrange = string(zminmax[0],format='(F4.2)')+','+string(zminmax[1],$
      format='(F4.2)')+','+strtrim(nzz,2)+','+strtrim(zlog,2)+' # [minz,maxz,dz,log?]'
    openw, lun, paramfile, /get_lun
    printf, lun, 'h100                 '+string(h100,format='(F4.2)')
    printf, lun, 'omega0               '+string(omega0,format='(F4.2)')
    printf, lun, 'omegal               '+string(omegal,format='(F4.2)')
    printf, lun, 'synthmodels          '+synthmodels
    printf, lun, 'imf                  '+imf
    printf, lun, 'sfhgrid              '+sfhgridstring
    printf, lun, 'redcurve             '+redcurvestring
    printf, lun, 'prefix               '+prefix
    printf, lun, 'redshift             '+zrange
    printf, lun, 'igm                  '+strtrim(igm,2)+' # [0=no, 1=yes]'
    printf, lun, 'maxold               0 # [0=no, 1=yes]'
    printf, lun, 'filterlist           '+strjoin(filters,',')
    free_lun, lun 

return
end

pro j0905_isedfit, supergrid, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber, debug=debug
; jm12feb08ucsd - derive stellar masses for J0905

    j0905path = j0905_path()
    isedpath = j0905_path(/isedfit)
    isedfit_sfhgrid_dir = j0905_path(/monte)
    sfhgrid_paramfile = getenv('HIZEA_DIR')+'/pro/j0905/j0905_sfhgrid.par'

    filters = j0905_filterlist()

    cat = mrdfits(j0905path+'j0905_photometry.fits.gz',1)
;   cat = mrdfits(j0905path+'j0905+5759_isedfit_input_v2.fits.gz',1)

; loop on each supergrid
    prefix = 'j0905'
    super = get_j0905_supergrid(supergrid,nsuper=nsuper)
    struct_print, super

    for gg = 0, nsuper-1 do begin
       splog, 'working on grid '+strtrim(super[gg].supergrid,2)

       paramfile = isedpath+prefix+'_supergrid'+string(super[gg].supergrid,$
         format='(i2.2)')+'_isedfit.par'
       j0905_write_paramfile, paramfile, prefix=prefix, super=super[gg]

; --------------------------------------------------
; build the models
       if keyword_set(models) then begin
          isedfit_models, paramfile, iopath=isedpath, clobber=clobber, $
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir
       endif

; --------------------------------------------------
; do the fitting!  
       if keyword_set(isedfit) then begin
          wise34 = where(strtrim(filters,2) eq 'wise_w3.par' or strtrim(filters,2) eq 'wise_w4.par')
          cat.ivarmaggies[wise34] = 0
          isedfit, paramfile, cat.maggies, cat.ivarmaggies, cat.z, iopath=isedpath, $
            clobber=clobber, sfhgrid_paramfile=sfhgrid_paramfile, $
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, outprefix=outprefix
       endif 

; --------------------------------------------------
; make a QAplot
       if keyword_set(qaplot) then begin
          yrange = [25,17]
          isedfit_qaplot, paramfile, result, iopath=isedpath, galaxy=cat.galaxy, $
            index=index, clobber=clobber, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, $
            outprefix=outprefix, xrange=xrange, yrange=yrange, /xlog
       endif
    endfor

return
end
