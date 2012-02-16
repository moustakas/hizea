pro build_j0905_sfhgrid, sfhgrid=sfhgrid, make_montegrid=make_montegrid, $
  clobber=clobber, debug=debug
; jm12feb08ucsd - build all the SFH grids we are going to need

    isedfit_sfhgrid_dir = j0905_path(/monte)
    sfhgrid_paramfile = getenv('HIZEA_DIR')+'/pro/j0905/j0905_sfhgrid.par'

; default grid
    build_isedfit_sfhgrid, '1', synthmodels='fsps', imf='chab', $
      redcurve=0, make_montegrid=make_montegrid, clobber=clobber, debug=debug, $
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, sfhgrid_paramfile=sfhgrid_paramfile
       
return
end
    
