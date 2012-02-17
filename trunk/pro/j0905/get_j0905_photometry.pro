pro get_j0905_photometry
; jm12feb16ucsd - gather all the photometry for J0905 and write out a
; FITS table

    path = j0905_path()
    
    filt = j0905_filterlist()
    weff = k_lambda_eff(filterlist=filt)
    nfilt = n_elements(filt)
    
; see http://cas.sdss.org/dr6/en/tools/explore/obj.asp?ra=136.34836&dec=57.986795
    phot = {galaxy: 'SDSS J090523.60+575912.4', ra: 136.3483648D, dec: 57.98679493D, $
      z: 0.7116, maggies: fltarr(nfilt), ivarmaggies: fltarr(nfilt)}
    
; galex; take the weighted mean of the GALEX NUV observations
    galex = im_read_tbl(path+'j0905_schiminovich_galex.txt',/silent)
    im_galex_to_maggies, galex, mm, ii, /psf

    nuv = im_weighted_mean(mm[1,*],error=1D/sqrt(ii[1,*]),wsigma=nuverr)
    phot.maggies[0] = reform(mm[0,3])
    phot.ivarmaggies[0] = reform(ii[0,3])
    phot.maggies[1] = nuv
    phot.ivarmaggies[1] = 1D/nuverr^2

; sdss
    sdss = mrdfits(hizea_path(/sdss)+'hizea_sdss_photo_dr7.fits.gz',1)
    spherematch, sdss.ra, sdss.dec, phot.ra, phot.dec, 1D/3600, m1, m2

    sdss_to_maggies, mm, ii, tsobj=sdss[m1], flux='model'
    phot.maggies[2:6] = mm
    phot.ivarmaggies[2:6] = ii

; irac
    minerr = [0.05,0.05]
    maggies = [144.0,113.3]*10^(-0.4*23.9)
    ivarmaggies = 1D/([1.7,1.2]*10^(-0.4*23.9))^2
    k_minerror, maggies, ivarmaggies, minerr

    phot.maggies[7:8] = maggies
    phot.ivarmaggies[7:8] = ivarmaggies

; synthesize photometry from the ground-based spectra
    lris = mrdfits(path+'spectra/j0905+5759_lris_flux_v120216.fits',1)
    amd = where(strmatch(filt,'*amd*'),namd,comp=other)
    maggies = k_project_filters(k_lambda_to_edges(lris.wavelength),$
      lris.flux,filterlist=filt[amd])
    maggieserr = k_project_filters(k_lambda_to_edges(lris.wavelength),$
      lris.error,filterlist=filt[amd])
    ivarmaggies = 1D/maggieserr^2
    phot.maggies[9:nfilt-1] = maggies
    phot.ivarmaggies[9:nfilt-1] = ivarmaggies
    
    djs_plot, lris.wavelength, 1D17*lris.flux, xsty=3, ysty=3
    djs_oplot, weff[amd], 1D17*phot.maggies[amd]*10^(-0.4*48.6)*im_light(/ang)/weff[amd]^2, $
      psym=6, symsize=4, color='green'
    djs_oplot, weff[other], 1D17*phot.maggies[other]*10^(-0.4*48.6)*im_light(/ang)/weff[other]^2, $
      psym=6, symsize=4, color='orange'
    
; write out
    minerr = replicate(0.05,nfilt)
    k_minerror, phot.maggies, phot.ivarmaggies, minerr
    
    im_mwrfits, phot, path+'j0905_photometry.fits', /clob
    
return
end
    
