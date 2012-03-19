pro get_sigmasfr_photometry, phot
; jm12mar19ucsd - gather all the photometry for SIGMASFR and write out a
; FITS table

    common sigmasfr, vagc
    
    path = sigmasfr_path()
    
    filt = sigmasfr_filterlist()
    weff = k_lambda_eff(filterlist=filt)
    nfilt = n_elements(filt)

    hst = rsex(path+'hst_sample.dat')
    ngal = n_elements(hst)
    
    phot = struct_trimtags(hst,select=['GALAXY','RA','DEC','Z'])
    phot = struct_addtags(phot,replicate({maggies: fltarr(nfilt), ivarmaggies: fltarr(nfilt)},ngal))
    
;; galex; take the weighted mean of the GALEX NUV observations
;    galex = im_read_tbl(path+'sigmasfr_schiminovich_galex.txt',/silent)
;    im_galex_to_maggies, galex, mm, ii, /psf
;
;    nuv = im_weighted_mean(mm[1,*],error=1D/sqrt(ii[1,*]),wsigma=nuverr)
;    phot.maggies[0] = reform(mm[0,3])
;    phot.ivarmaggies[0] = reform(ii[0,3])
;    phot.maggies[1] = nuv
;    phot.ivarmaggies[1] = 1D/nuverr^2

; sdss
    if (n_elements(vagc) eq 0L) then vagc = $
      mrdfits(getenv('VAGC_REDUX')+'/object_sdss_imaging.fits.gz',1)
    spherematch, vagc.ra, vagc.dec, phot.ra, phot.dec, 1D/3600, m1, m2
    sdss_to_maggies, mm, ii, calib=vagc[m1], flux='model'

    phot[m2].maggies[2:6] = mm
    phot[m2].ivarmaggies[2:6] = ii

; wise
    wise = mrdfits(path+'hst_wise_allsky.fits.gz',1)
    spherematch, wise.ra, wise.dec, phot.ra, phot.dec, 5D/3600, m1, m2
    wise_to_maggies, wise[m1], mm, ii, /mpro
    
    phot[m2].maggies[7:10] = mm
    phot[m2].ivarmaggies[7:10] = ii

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
    
