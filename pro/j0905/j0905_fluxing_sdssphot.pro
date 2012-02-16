;pro j0905_fluxing_sdssphot
;+
;  amd111213 -- written to scale Keck/LRIS spectrum for J0905 to SDSS
;               photometry
;
;-

keck = mrdfits('j0905+5759_lris_110307.fits', 1)

mmt_a = rd1dspec('T12a_J0905+5759_006.ms.fits')
mmt_b = rd1dspec('T12b_J0905+5759_006.ms.fits')

sdss_flux1 = mrdfits('spSpec-51902-0483-601.fit', 0, hdr)
sdss_wave1 = 10.0^(findgen(sxpar(hdr, 'NAXIS1')) * sxpar(hdr, 'COEFF1') + $
            sxpar(hdr, 'COEFF0')) 

sdss_flux2 = mrdfits('spSpec-51924-0483-628.fit', 0, hdr)
sdss_wave2 = 10.0^(findgen(sxpar(hdr, 'NAXIS1')) * sxpar(hdr, 'COEFF1') + $
            sxpar(hdr, 'COEFF0')) 


sdss_flux3 = mrdfits('spSpec-51942-0483-639.fit', 0, hdr)
sdss_wave3 = 10.0^(findgen(sxpar(hdr, 'NAXIS1')) * sxpar(hdr, 'COEFF1') + $
            sxpar(hdr, 'COEFF0')) 

;-----------------------------------------------------------------------------

model_mag = [19.912359, 19.560743, 19.384109, 19.071833, 19.065737]
psf_mag = [19.919191, 19.569157, 19.391365, 19.089397, 19.073381]
psf_err = [0.042183, 0.028049, 0.021074, 0.021189, 0.073866]

ugriz_lam0 = [3540., 4770., 6222., 7632., 9049]
ugriz_dlam = [463, 988, 955, 1064, 1248]

photo_flux = ugriz_mag2flam(psf_mag, psf_err, flux_err = photo_err)

;-------------------------------------------------------------------------

kscale = fltarr(n_elements(ugriz_lam0))
kerr = fltarr(n_elements(ugriz_lam0))
for i=0L, n_elements(ugriz_lam0)-1L do begin
  rel = where(keck.wavelength lt ugriz_lam0[i]+ugriz_dlam[i]/2. and $
          keck.wavelength gt ugriz_lam0[i]-ugriz_dlam[i]/2. and $
          (keck.wavelength lt 7550. or keck.wavelength gt 7700.), nrel)  
  kflux = amd_wav(keck[rel].flux, keck[rel].error) 
  kscale[i] = photo_flux[i]*1.d17 / kflux.value
  kerr[i] = kscale[i]*sqrt((photo_err[i]/photo_flux[i])^2+$
    (kflux.sigma/kflux.value)^2)
endfor

splog, kscale, kerr
;kerr[1] = kerr[1]*5.

dfpsplot, 'j0905_fluxing_sdssphot.ps', /landscape, /color


  plot, ugriz_lam0, kscale, psym=6
  djs_oploterr, ugriz_lam0, kscale, xerr = ugriz_dlam/2, $
    yerr=kerr, psym=1, sym=2, thick=6

; fit a polynomial to the above
  poly_expr = 'P(0)*x + P(1)'
  poly_p0 = [1., 1.]
  ;poly_expr = 'P(0)*x^2 + P(1)*x + P(2)'
  ;poly_p0 = [1., 1., 1.]
  poly_p = mpfitexpr(poly_expr, ugriz_lam0, kscale, kerr, poly_p0, /quiet)
  poly_fit = mpevalexpr(poly_expr, ugriz_lam0, poly_p)
  oplot, ugriz_lam0, poly_fit, color=djs_icolor('red')  
  scale = mpevalexpr(poly_expr, keck.wavelength, poly_p)

  replace = where(scale gt 1.3, nreplace)
  scale[replace]=1.3
  replace = where(scale lt 1.05, nreplace)
  scale[replace]=1.05

  oplot, keck.wavelength, scale

  plot, keck.wavelength, keck.flux*scale, $
    /xs, yr=[0, 12], xr=[4100, 8900], $
    xtitle = 'Observered Wavelength', charthick=2
  oplot, keck.wavelength, keck.flux, color=djs_icolor('red')
  djs_oploterr, ugriz_lam0, photo_flux*1e17, xerr = ugriz_dlam/2, $
      yerr=photo_err*1e17, psym=1, sym=2, thick=6

  plot, keck.wavelength, keck.flux*scale, $
    /xs, yr=[0, 12], xr=[4100, 8900], $
    xtitle = 'Observered Wavelength', charthick=2
  oplot, keck.wavelength, keck.flux, color=djs_icolor('red')
  djs_oploterr, ugriz_lam0, photo_flux*1e17, xerr = ugriz_dlam/2, $
    yerr=photo_err*1e17, psym=1, sym=2, thick=6

  djs_oplot, sdss_wave1, medsmooth(sdss_flux1[*,0],9), color = 'red'
  djs_oplot, sdss_wave2, medsmooth(sdss_flux2[*,0],9), color = 'orange'
  djs_oplot, sdss_wave3, medsmooth(sdss_flux3[*,0],9), color = 'cyan'

  djs_oplot, mmt_a.wave, medsmooth(mmt_a.spec*1e17 /1.8, 9), color='blue'
  djs_oplot, mmt_b.wave, medsmooth(mmt_b.spec*1e17 /2.2, 9), color='green'

dfpsclose

stop



end
