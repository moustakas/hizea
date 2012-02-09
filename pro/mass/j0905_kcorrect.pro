pro j0905_kcorrect, noirac=noirac, clobber=clobber
; jm12feb08ucsd - compute K-corrections for J0905

;   kcorr = mrdfits(j0905path+'j0905_kcorrect.fits.gz',1)
;   kcorrect_qaplot, kcorr, psfile=j0905path+'qaplot_j0905_kcorrect.ps', $
;     in_filterlist=j0905_filterlist(), /clobber, vname='default.nolines'
    
    j0905path = j0905_path()
    h100 = 0.7
    vname = 'default.nolines'

    if keyword_set(noirac) then suffix = '_noirac' else suffix = ''
    kcorrfile = j0905path+'j0905_kcorrect'+suffix+'.fits'
    if im_file_test(kcorrfile+'.gz',clobber=clobber) then return

    cat = mrdfits(j0905path+'j0905+5759_isedfit_input.fits.gz',1)
    maggies = cat.maggies
    ivarmaggies = cat.ivarmaggies
    zobj = cat.z
    filterlist = j0905_filterlist()
    if keyword_set(noirac) then begin
       toss = where(strmatch(filterlist,'*irac*'))
       ivarmaggies[toss,*] = 0.0
       outprefix = 'j0905_noirac'
    endif
    ngal = n_elements(zobj)
    nfilt = n_elements(filterlist)

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
      k_fnuvugrizjhk_absmag_01:         fltarr(10)-999.0,$
      k_fnuvugrizjhk_absmag_ivar_01:    fltarr(10)-999.0,$
      k_fnuvugrizjhk_kcorrect_01:       fltarr(10)-999.0}
    kcorr = replicate(kcorr,ngal)

    kcorr.k_zobj = zobj
    kcorr.k_maggies = maggies
    kcorr.k_ivarmaggies = ivarmaggies
    
; compute k-corrections
    out_filterlist = [galex_filterlist(),sdss_filterlist(),twomass_filterlist()]
    
    splog, 'Computing [FN]UV-ugriz-JHK K-corrections'
    kcorrect_01 = im_kcorrect(zobj,maggies,ivarmaggies,$
      filterlist,out_filterlist,band_shift=0.1,chi2=chi2,mass=mass,$
      coeffs=coeffs,bestmaggies=bestmaggies,absmag=absmag_01,$
      ivarabsmag=absmag_ivar_01,uvflux=uvflux,$;clineflux=cflux,uvflux=uvflux,$
      /silent,vname=vname,h100=h100) ; AB, band_shift=0.1

; get the apparent magnitudes in ugrizJHKs
    k_load_vmatrix, vmatrix, lambda, vname=vname
    k_reconstruct_maggies, coeffs, zobj, appmaggies, vmatrix=vmatrix, $
      lambda=lambda, filterlist=out_filterlist
    kcorr.k_mobs_fnuvugrizjhk = reform(-2.5*alog10(appmaggies))
    
    kcorr.k_bestmaggies = bestmaggies
    kcorr.k_mass = alog10(mass) ; h=0.7, Chabrier
    kcorr.k_coeffs = coeffs
    kcorr.k_chi2 = chi2
    kcorr.k_uvflux = uvflux
    kcorr.k_fnuvugrizjhk_absmag_01      = absmag_01
    kcorr.k_fnuvugrizjhk_absmag_ivar_01 = absmag_ivar_01
    kcorr.k_fnuvugrizjhk_kcorrect_01    = kcorrect_01

    im_mwrfits, kcorr, kcorrfile, /clobber

return
end
