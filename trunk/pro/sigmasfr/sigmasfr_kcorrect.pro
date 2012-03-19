pro sigmasfr_kcorrect, clobber=clobber
; jm12mar19ucsd - compute K-corrections for the SIGMASFR project

;   kcorr = mrdfits(sigmasfrpath+'sigmasfr_kcorrect.fits.gz',1)
;   kcorrect_qaplot, kcorr, psfile=sigmasfrpath+'qaplot_sigmasfr_kcorrect.ps', $
;     in_filterlist=sigmasfr_filterlist(), /clobber, vname='default.nolines'
    
    sigmasfrpath = sigmasfr_path()

    h100 = 0.7
    vname = 'default'
;   vname = 'default.nolines'

    kcorrfile = sigmasfrpath+'sigmasfr_kcorrect.fits'
    if im_file_test(kcorrfile+'.gz',clobber=clobber) then return

    filters = sigmasfr_filterlist()
    nfilt = n_elements(filters)

    cat = mrdfits(sigmasfrpath+'sigmasfr_photometry.fits.gz',1)
    ngal = n_elements(cat)

    wise34 = where(strtrim(filters,2) eq 'wise_w3.par' or strtrim(filters,2) eq 'wise_w4.par')
    ivarmaggies = cat.ivarmaggies
    ivarmaggies[wise34] = 0

    kcorr = {$
      k_zobj:                          -999.0,$
      k_maggies:                fltarr(nfilt),$
      k_ivarmaggies:            fltarr(nfilt),$
      k_bestmaggies:            fltarr(nfilt),$
      k_mass:                          -999.0,$
      k_coeffs:                     fltarr(5),$
      k_chi2:                          -999.0,$
      k_uvflux:                     fltarr(2),$

      k_mobs_fnuvugrizjhk:              fltarr(10)-999.0,$
      k_fnuvugrizjhk_absmag_00:         fltarr(10)-999.0,$
      k_fnuvugrizjhk_absmag_ivar_00:    fltarr(10)-999.0,$
      k_fnuvugrizjhk_kcorrect_00:       fltarr(10)-999.0}
    kcorr = replicate(kcorr,ngal)

    kcorr.k_zobj = cat.z
    kcorr.k_maggies = cat.maggies
    kcorr.k_ivarmaggies = cat.ivarmaggies
    
; compute k-corrections
    out_filters = [galex_filterlist(),sdss_filterlist(),twomass_filterlist()]
    
    splog, 'Computing [FN]UV-ugriz-JHK K-corrections'
    kcorrect_00 = im_kcorrect(cat.z,cat.maggies,ivarmaggies,$
      filters,out_filters,band_shift=0.0,chi2=chi2,mass=mass,$
      coeffs=coeffs,bestmaggies=bestmaggies,absmag=absmag_00,$
      ivarabsmag=absmag_ivar_00,uvflux=uvflux,$;clineflux=cflux,$
      /silent,vname=vname,h100=h100) ; AB, band_shift=0.0

; get the apparent and absolute magnitudes in ugrizJHKs
    k_load_vmatrix, vmatrix, lambda, vname=vname
    k_reconstruct_maggies, coeffs, cat.z, appmaggies, vmatrix=vmatrix, $
      lambda=lambda, filterlist=out_filters
    kcorr.k_mobs_fnuvugrizjhk = reform(-2.5*alog10(appmaggies))
    
    kcorr.k_bestmaggies = bestmaggies
    kcorr.k_mass = alog10(mass) ; h=0.7, Chabrier
    kcorr.k_coeffs = coeffs
    kcorr.k_chi2 = chi2
    kcorr.k_uvflux = uvflux
    kcorr.k_fnuvugrizjhk_absmag_00      = absmag_00
    kcorr.k_fnuvugrizjhk_absmag_ivar_00 = absmag_ivar_00
    kcorr.k_fnuvugrizjhk_kcorrect_00    = kcorrect_00

    im_mwrfits, kcorr, kcorrfile, /clobber

return
end
