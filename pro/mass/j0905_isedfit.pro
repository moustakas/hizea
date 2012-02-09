pro j0905_write_paramfile, paramfile, sfhgrid=sfhgrid, prefix=prefix, $
  zminmax=zminmax, nzz=nzz, filters=filters, igm=igm

    h100 = 0.7
    omega0 = 0.3
    omegal = 0.7

    synthmodels = 'bc03'
    imf = 'chab'
    zlog = '0'
    redcurvestring = redcurve2string(1)
    sfhgridstring = strtrim(sfhgrid,2)
    
    splog, 'Writing '+paramfile
    zrange = string(zminmax[0],format='(F4.2)')+','+string(zminmax[1],$
      format='(F4.2)')+','+strtrim(nzz,2)+','+zlog+' # [minz,maxz,dz,log?]'
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

pro j0905_isedfit, models=models, write_paramfile=write_paramfile, $
  isedfit=isedfit, qaplot=qaplot, clobber=clobber, debug=debug, $
  noirac=noirac, j0905=j0905
; jm12feb08ucsd - derive stellar masses for J0905

    j0905path = j0905_path()
    isedpath = j0905_path(/isedfit)
    isedfit_sfhgrid_dir = j0905_path(/monte)
    sfhgrid_paramfile = getenv('HIZEA_DIR')+'/pro/mass/j0905_sfhgrid.par'

    prefix = 'j0905_sfhgrid01'
    paramfile = isedpath+prefix+'_isedfit.par'
    filters = j0905_filterlist()
    zminmax = [0.71,0.72]
    nzz = 3
    igm = 1
    sfhgrid = 1
    
    j0905_write_paramfile, paramfile, prefix=prefix, sfhgrid=sfhgrid, $
      zminmax=zminmax, nzz=nzz, filters=filters, igm=igm

; --------------------------------------------------
; build the models
    if keyword_set(models) then begin
       isedfit_models, paramfile, iopath=isedpath, clobber=clobber, $
         isedfit_sfhgrid_dir=isedfit_sfhgrid_dir
    endif

; --------------------------------------------------
; do the fitting!  
    if keyword_set(isedfit) then begin
       cat = mrdfits(j0905path+'j0905+5759_isedfit_input.fits.gz',1)

       isedfit, paramfile, cat.maggies, cat.ivarmaggies, cat.z, iopath=isedpath, $
         clobber=clobber, sfhgrid_paramfile=sfhgrid_paramfile, $
         isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, outprefix=outprefix
    endif 

; --------------------------------------------------
; make a QAplot
    if keyword_set(qaplot) then begin
       cat = mrdfits(j0905path+'j0905+5759_isedfit_input.fits.gz',1)
       galaxy = 'J0905+5759'
;      if keyword_set(noirac) then outprefix = 'hizea_noirac'

       yrange = [25,17]
       isedfit_qaplot, paramfile, result, iopath=isedpath, galaxy=galaxy, $
         index=index, clobber=clobber, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, $
         outprefix=outprefix, xrange=xrange, yrange=yrange, /xlog
    endif

return
end
