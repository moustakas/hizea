pro build_j0905_supergrid, supergrid, make_montegrid=make_montegrid, clobber=clobber
; jm12feb16ucsd - build all the SFH grids we are going to need using
; the supergrid parameter file

; echo "build_j0905_supergrid, [1,2], /make, /cl" | idl > & ~/j0905.1.log &     
; echo "build_j0905_supergrid, [3,4], /make, /cl" | idl > & ~/j0905.2.log &

    isedfit_sfhgrid_dir = j0905_path(/monte)
    sfhgrid_paramfile = getenv('HIZEA_DIR')+'/pro/j0905/j0905_sfhgrid.par'

    super = get_j0905_supergrid(supergrid,nsuper=nsuper)
    struct_print, super

    for ii = 0, nsuper-1 do begin
       build_isedfit_sfhgrid, super[ii].sfhgrid, synthmodels=strtrim(super[ii].synthmodels,2), $
         imf=strtrim(super[ii].imf,2), redcurve=super[ii].redcurve, make_montegrid=make_montegrid, $
         clobber=clobber, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, sfhgrid_paramfile=sfhgrid_paramfile
    endfor

return
end
