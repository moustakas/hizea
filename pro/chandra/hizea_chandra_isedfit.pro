pro hizea_chandra_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  kcorrect=kcorrect, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, clobber=clobber
; jm13sep01siena - fit Paul's Chandra sample
    
;   prefix = 'test'
    prefix = 'hizea_chandra'
    isedfit_dir = getenv('HIZEA_DATA')+'/chandra/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

    filterlist = sigmasfr_filterlist()
    
; --------------------------------------------------
; write the parameter file
    if keyword_set(write_paramfile) then begin
       cat = mrdfits(isedfit_dir+'hizea_chandra_photometry.fits.gz',1)
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, spsmodels='fsps_v2.4_miles', $
         imf='chab', redcurve='charlot', /igm, use_redshift=cat.z, $
         nmodel=20000L, age=[1.0,9.0], tau=[0.1,10.0], Zmetal=[0.005,0.03], $
         AV=[0.0,5.0], mu=[0.0,0.7], /flatAV, /flatmu, pburst=0.5, $
         interval_pburst=8.0, tburst=[1.0,9.0], fburst=[0.1,5.0], $
         dtburst=[0.01,1.0], trunctau=[0.005,0.2], fractrunc=0.8, $
         oiiihb=[-1.0,1.0], /flatfburst, /flatdtburst, /nebular, $
         clobber=clobber
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, spsmodels='fsps_v2.4_miles', $
         imf='chab', redcurve='smc', /igm, use_redshift=cat.z, $
         nmodel=20000L, age=[1.0,9.0], tau=[0.1,10.0], Zmetal=[0.005,0.03], $
         AV=[0.0,5.0], mu=[0.0,0.7], /flatAV, /flatmu, pburst=0.5, $
         interval_pburst=8.0, tburst=[1.0,9.0], fburst=[0.1,5.0], $
         dtburst=[0.01,1.0], trunctau=[0.005,0.2], fractrunc=0.8, $
         oiiihb=[-1.0,1.0], /flatfburst, /flatdtburst, /nebular, $
         clobber=clobber, /append
    endif

; --------------------------------------------------
; build the Monte Carlo grids    
    if keyword_set(build_grids) then begin
       isedfit_montegrids, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, clobber=clobber
    endif

; --------------------------------------------------
; calculate the model photometry 
    if keyword_set(model_photometry) then begin
       isedfit_models, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber
    endif

; --------------------------------------------------
; generate the model photometry QAplots
    if keyword_set(qaplot_models) then begin
       cat = mrdfits(isedfit_dir+'hizea_chandra_photometry.fits.gz',1)
       thesefilters = ['galex_NUV','sdss_g0','sdss_r0','sdss_i0','wise_w1']
       isedfit_qaplot_models, isedfit_paramfile, cat.maggies, $
         cat.ivarmaggies, cat.z, isedfit_dir=isedfit_dir, $
         thissfhgrid=thissfhgrid, thesefilters=thesefilters, clobber=clobber
    endif
    
; --------------------------------------------------
; fit! drop the [ch3,4] WISE photometry
    if keyword_set(isedfit) then begin
       cat = mrdfits(isedfit_dir+'hizea_chandra_photometry.fits.gz',1)
       wise34 = where(strtrim(filterlist,2) eq 'wise_w3.par' or $
         strtrim(filterlist,2) eq 'wise_w4.par')
       cat.ivarmaggies[wise34] = 0       
       isedfit, isedfit_paramfile, cat.maggies, cat.ivarmaggies, $
         cat.z, ra=cat.ra, dec=cat.dec, isedfit_dir=isedfit_dir, $
         thissfhgrid=thissfhgrid, clobber=clobber
    endif 

; --------------------------------------------------
; compute K-corrections
    if keyword_set(kcorrect) then begin
       isedfit_kcorrect, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         absmag_filterlist=bessell_filterlist(), band_shift=0.0, $
         clobber=clobber
    endif 

; --------------------------------------------------
; generate spectral energy distribution (SED) QAplots
    if keyword_set(qaplot_sed) then begin
       isedfit_qaplot_sed, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber, /xlog, xrange=[1D3,8D4], yrange=[26,13]
    endif

return
end
