pro get_sigmasfr_photometry, phot
; jm12mar19ucsd - gather all the photometry for SIGMASFR and write out a
; FITS table
;
; ad12mar19ucsd - added in pre-schiminovich galex photometry
;

    common sigmasfr, vagc
    
    path = sigmasfr_path()
    
    filt = sigmasfr_filterlist()
    weff = k_lambda_eff(filterlist=filt)
    nfilt = n_elements(filt)

    hst = rsex(path+'hst_sample.dat')
    ngal = n_elements(hst)
    
    phot = struct_trimtags(hst,select=['GALAXY','RA','DEC','Z'])
    phot = struct_addtags(phot,replicate({maggies: fltarr(nfilt), $
      ivarmaggies: fltarr(nfilt)},ngal))
    
;; galex; take the weighted mean of the GALEX NUV observations
;    galex = im_read_tbl(path+'sigmasfr_schiminovich_galex.txt',/silent)
;    im_galex_to_maggies, galex, mm, ii, /psf
;
;    nuv = im_weighted_mean(mm[1,*],error=1D/sqrt(ii[1,*]),wsigma=nuverr)
;    phot.maggies[0] = reform(mm[0,3])
;    phot.ivarmaggies[0] = reform(ii[0,3])
;    phot.maggies[1] = nuv
;    phot.ivarmaggies[1] = 1D/nuverr^2

; galex (pre-schiminovich)
    galex = im_read_tbl(path+'galexview_output.tbl',/silent)
    im_galex_to_maggies, galex, mm, ii, /psf
; loop through the GALEX table
    for jj=0L, n_elements(phot)-1L do begin
       match=where(galex.UPLOADSEARCHID eq jj+1L, nmatch)
; find useful FUV data
       fuvgood=where(mm[0,match] gt 0., nfuvgood)
       if nfuvgood lt 1L then begin
          phot[jj].maggies[0] = 0.
          phot[jj].ivarmaggies[0] = 0.
       endif else begin
          if nfuvgood eq 1L then begin
             phot[jj].maggies[0] = mm[0,match[fuvgood]]
             phot[jj].ivarmaggies[0] = ii[0,match[fuvgood]]      
          endif else begin
             fuv = im_weighted_mean(mm[0,match[fuvgood]],$
               error=1D/sqrt(ii[0,match[fuvgood]]),wsigma=fuverr)
             phot[jj].maggies[0] = fuv
             phot[jj].ivarmaggies[0] = 1D/fuverr^2
          endelse
       endelse
; find useful NUV data
       nuvgood=where(mm[1,match] gt 0., nnuvgood)
       if nnuvgood lt 1L then begin
          phot[jj].maggies[1] = 0.
          phot[jj].ivarmaggies[1] = 0.
       endif else begin
          if nnuvgood eq 1L then begin
             phot[jj].maggies[1] = mm[1,match[nuvgood]]
             phot[jj].ivarmaggies[1] = ii[1,match[nuvgood]]      
          endif else begin
             nuv = im_weighted_mean(mm[1,match[nuvgood]],$
               error=1D/sqrt(ii[1,match[nuvgood]]),wsigma=nuverr)
             phot[jj].maggies[1] = nuv
             phot[jj].ivarmaggies[1] = 1D/nuverr^2
          endelse
       endelse  
    endfor

; sdss
    if (n_elements(vagc) eq 0L) then vagc = $
      mrdfits(getenv('VAGC_REDUX')+'/object_sdss_imaging.fits.gz',1)
    spherematch, vagc.ra, vagc.dec, phot.ra, phot.dec, 1D/3600, m1, m2
    sdss_to_maggies, mm, ii, calib=vagc[m1], flux='model'

    phot[m2].maggies[2:6] = mm
    phot[m2].ivarmaggies[2:6] = ii

; spitzer
    spit = rsex(path+'hst_sample_spitzer.dat')
    match, strtrim(spit.galaxy,2), strtrim(phot.galaxy,2), m1, m2
    ww = where(spit[m1].ch1_flux gt 0)
    phot[m2[ww]].maggies[7:8] = transpose([[spit[m1[ww]].ch1_flux],[spit[m1[ww]].ch2_flux]]*10^(-0.4*23.9))
    phot[m2[ww]].ivarmaggies[7:8] = 1D/(transpose([[spit[m1[ww]].ch1_sig],[spit[m1[ww]].ch2_sig]]*10^(-0.4*23.9)))^2

; wise
    wise = mrdfits(path+'hst_wise_allsky.fits.gz',1)
    spherematch, wise.ra, wise.dec, phot.ra, phot.dec, 5D/3600, m1, m2
    wise_to_maggies, wise[m1], mm, ii, /mpro
    
    phot[m2].maggies[9:12] = mm
    phot[m2].ivarmaggies[9:12] = ii

; write out
;   minerr = replicate(0.05,nfilt)
;   maggies = phot.maggies
;   ivarmaggies = phot.ivarmaggies
;   k_minerror, maggies, ivarmaggies, minerr
;   phot.maggies = maggies
;   phot.ivarmaggies = ivarmaggies
    im_mwrfits, phot, path+'sigmasfr_photometry.fits', /clob
    
return
end
    
