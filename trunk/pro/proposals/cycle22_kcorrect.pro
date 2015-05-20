pro cycle22_kcorrect
; jm14apr04siena - compute K-corrections for our Cycle 22 HST proposal

    gal = [$
      'J1506+5402',$
      'J0905+5759',$
      'J1341-0321',$
      'J0944+0930',$
      'J2140+1209',$
      'J0826+4305',$
      'J1613+2834',$
      'J1219+0336',$
      'J1107+0417',$
      'J0901+0314',$
      'J0106-1023',$
      'J2116-0634']
    ngal = n_elements(gal)
    
    sigmasfrpath = sigmasfr_path()
    path = getenv('HIZEA_DATA')+'/'

    h100 = 0.7
    vname = 'default'
;   vname = 'default.nolines'

    kcorrfile = path+'cycle22_kcorrect.fits'
    if im_file_test(kcorrfile+'.gz',clobber=clobber) then return

    filters = sigmasfr_filterlist()
    nfilt = n_elements(filters)

    cat = mrdfits(sigmasfrpath+'sigmasfr_photometry.fits.gz',1)
    cat[10].galaxy = 'J0106-1023' ; typo!

    match, gal, strtrim(cat.galaxy,2), m1, m2
    srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
    cat = cat[m2]

    wise34 = where(strtrim(filters,2) eq 'wise_w3.par' or $
      strtrim(filters,2) eq 'wise_w4.par')
    ivarmaggies = cat.ivarmaggies
    ivarmaggies[wise34] = 0

    kcorr = {$
      galaxy:                            '',$
      z:                          -999.0,$
      maggies:                fltarr(nfilt),$
      ivarmaggies:            fltarr(nfilt),$
      bestmaggies:            fltarr(nfilt),$
      mass:                          -999.0,$
      coeffs:                     fltarr(5),$
      chi2:                          -999.0,$
      uvflux:                     fltarr(2),$

      uvj_absmag_00:         fltarr(3)-999.0,$
      uvj_absmag_ivar_00:    fltarr(3)-999.0,$
      uvj_kcorrect_00:       fltarr(3)-999.0}
    kcorr = replicate(kcorr,ngal)

    kcorr.galaxy = cat.galaxy
    kcorr.z = cat.z
    kcorr.maggies = cat.maggies
    kcorr.ivarmaggies = cat.ivarmaggies
    
; compute k-corrections
    out_filters = ['bessell_U.par','bessell_V.par','twomass_J.par']
    
    splog, 'Computing K-corrections'
    kcorrect_00 = im_kcorrect(cat.z,cat.maggies,ivarmaggies,$
      filters,out_filters,band_shift=0.0,chi2=chi2,mass=mass,$
      coeffs=coeffs,bestmaggies=bestmaggies,absmag=absmag_00,$
      ivarabsmag=absmag_ivar_00,uvflux=uvflux,$;clineflux=cflux,$
      /silent,vname=vname,h100=h100) ; AB, band_shift=0.0

; get the apparent and absolute magnitudes in ugrizJHKs
    k_load_vmatrix, vmatrix, lambda, vname=vname
    k_reconstruct_maggies, coeffs, cat.z, appmaggies, vmatrix=vmatrix, $
      lambda=lambda, filterlist=out_filters
;   kcorr.k_mobs_uvj = reform(-2.5*alog10(appmaggies))
    
    kcorr.bestmaggies = bestmaggies
    kcorr.mass = alog10(mass) ; h=0.7, Chabrier
    kcorr.coeffs = coeffs
    kcorr.chi2 = chi2
    kcorr.uvflux = uvflux
    kcorr.uvj_absmag_00      = absmag_00
    kcorr.uvj_absmag_ivar_00 = absmag_ivar_00
    kcorr.uvj_kcorrect_00    = kcorrect_00

    im_mwrfits, kcorr, kcorrfile, /clobber

return
end
    
