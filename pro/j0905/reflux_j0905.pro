function read_sdssspec, file
    flux = mrdfits(file,0,hdr)
    wave = 10.0^(findgen(sxpar(hdr,'NAXIS1'))*sxpar(hdr, 'COEFF1') + $
      sxpar(hdr,'COEFF0')) 
    sdss = {wave: wave, flux: reform(flux[*,0])*1D-17, $
      ferr: reform(flux[*,1])*1D-17}
return, sdss
end

pro reflux_j0905
; jm12feb16ucsd - adjust the fluxing of the LRIS spectrum of J0905

    path = j0905_path()+'spectra/'
    ra = 136.3483648D
    dec = 57.98679493D
  
    glactc, ra, dec, 2000.0, gl, gb, 1, /deg
    ebv = dust_getval(gl,gb,/interp,/noloop)

; read the SDSS photometry    
    phot = mrdfits(hizea_path(/sdss)+'hizea_sdss_photo_dr7.fits.gz',1)
    spherematch, phot.ra, phot.dec, ra, dec, 1D/3600, m1, m2

    sdss_to_maggies, ugriz, tsobj=phot[m1], flux='model'
    weff = k_lambda_eff()
    ugriz_flam = ugriz*10^(-0.4*48.6)*im_light(/ang)/weff^2
    
    normindx = 3 ; =iband
    normfilt = 'sdss_i0.par'
    
; Keck/LRIS spectrum
    keck = mrdfits(path+'j0905+5759_lris_110307.fits', 1)
;   keck1 = mrdfits(path+'j0905+5759_lris_110307.fits', 1)
;   keep = where(keck1.wavelength gt 3950.0,nkeep)
;   keck = im_empty_structure(keck1,ncopies=nkeep)
;   keck.wavelength = keck1[keep].wavelength
;   keck.flux = keck1[keep].flux
;   keck.error = keck1[keep].error

    dust = 10^(0.4*ebv*k_lambda(keck.wavelength,/odon))
    keck.flux = keck.flux*1D-17*dust
    keck.error = keck.error*1D-17*dust
    
    mag = k_project_filters(k_lambda_to_edges(keck.wavelength),keck.flux,$
      filterlist=normfilt)
    keck.flux = keck.flux*(ugriz[normindx]/mag[0])
    keck.error = keck.error*(ugriz[normindx]/mag[0])
    newkeck = keck
    
    keck_ugriz = k_project_filters(k_lambda_to_edges(keck.wavelength),keck.flux)
    keck_ugriz_flam = keck_ugriz*10^(-0.4*48.6)*im_light(/ang)/weff^2
    
; MMT spectra
    mmt_a = rd1dspec('T12a_J0905+5759_006.ms.fits',datapath=path)
    mag = k_project_filters(k_lambda_to_edges(mmt_a.wave),mmt_a.spec,$
      filterlist=normfilt)
    mmt_a.spec = mmt_a.spec*(ugriz[normindx]/mag[0])
    mmt_a.sigspec = mmt_a.sigspec*(ugriz[normindx]/mag[0])

    mmt_ugriz = k_project_filters(k_lambda_to_edges(mmt_a.wave),mmt_a.spec)
    mmt_ugriz_flam = mmt_ugriz*10^(-0.4*48.6)*im_light(/ang)/weff^2
    
    mmt_b = rd1dspec('T12b_J0905+5759_006.ms.fits',datapath=path)
    mag = k_project_filters(k_lambda_to_edges(mmt_b.wave),mmt_b.spec,$
      filterlist=normfilt)
    mmt_b.spec = mmt_b.spec*(ugriz[normindx]/mag[0])
    mmt_b.sigspec = mmt_b.sigspec*(ugriz[normindx]/mag[0])

; SDSS spectra    
    sdss1 = read_sdssspec(path+'spSpec-51902-0483-601.fit')
    mag = k_project_filters(k_lambda_to_edges(sdss1.wave),sdss1.flux,filterlist=normfilt)
    sdss1.flux = sdss1.flux*(ugriz[normindx]/mag[0])
    sdss1.ferr = sdss1.ferr*(ugriz[normindx]/mag[0])

    sdss2 = read_sdssspec(path+'spSpec-51924-0483-628.fit')
    mag = k_project_filters(k_lambda_to_edges(sdss2.wave),sdss2.flux,filterlist=normfilt)
    sdss2.flux = sdss2.flux*(ugriz[normindx]/mag[0])
    sdss2.ferr = sdss2.ferr*(ugriz[normindx]/mag[0])

    sdss3 = read_sdssspec(path+'spSpec-51942-0483-639.fit')
    mag = k_project_filters(k_lambda_to_edges(sdss3.wave),sdss3.flux,filterlist=normfilt)
    sdss3.flux = sdss3.flux*(ugriz[normindx]/mag[0])
    sdss3.ferr = sdss3.ferr*(ugriz[normindx]/mag[0])

; fit a polynomial to the ugriz photometry
    sdsscoeff = poly_fit(weff[0:2],ugriz_flam[0:2],2)
    keckcoeff = poly_fit(weff[0:2],keck_ugriz_flam[0:2],2)
    keckcorr = (poly(keck.wavelength,sdsscoeff)/poly(keck.wavelength,keckcoeff))>1
    newkeck.flux = keck.flux*keckcorr

    newkeck_ugriz = k_project_filters(k_lambda_to_edges(newkeck.wavelength),newkeck.flux)
    newkeck_ugriz_flam = newkeck_ugriz*10^(-0.4*48.6)*im_light(/ang)/weff^2

    niceprint, weff, -2.5*alog10(newkeck_ugriz/ugriz), -2.5*alog10(keck_ugriz/ugriz)

;   mmtc = im_fitcontinuum(mmt_a.wave,mmt_a.spec,mmt_a.sigspec,$
;     method=2,/tell,medwidth=201,smoothwidth=201)
;   keckc = im_fitcontinuum(keck.wavelength,keck.flux,keck.error,$
;     method=2,/tell,medwidth=201,smoothwidth=201)
    
    im_mwrfits, newkeck, path+'j0905+5759_lris_flux_v120216.fits', /clobber, /nogzip
    
; make a QAplot
    psfile = path+'qa_reflux_j0905.ps'
    im_plotconfig, 0, pos, psfile=psfile

    djs_plot, [0], [0], /nodata, yrange=[0,15], xsty=3, $
      ysty=3, xrange=[2900,10000], position=pos
    djs_oplot, keck.wavelength, medsmooth(keck.flux*1D17,10), psym=10, color='green'
    djs_oplot, newkeck.wavelength, medsmooth(newkeck.flux*1D17,10), psym=10, color='red'

;   djs_oplot, sdss1.wave, medsmooth(sdss1.flux*1D17,14), psym=10, color='red'
;   djs_oplot, sdss2.wave, medsmooth(sdss2.flux*1D17,14), psym=10, color='blue'
;   djs_oplot, sdss3.wave, medsmooth(sdss3.flux*1D17,14), psym=10, color='purple'
;   djs_oplot, mmt_a.wave, medsmooth(mmt_a.spec*1D17,10), psym=10, color='orange'
;   djs_oplot, mmt_b.wave, medsmooth(mmt_b.spec*1D17,10), psym=10, color='cyan'

;   djs_oplot, weff, 1D17*mmt_ugriz_flam, psym=7, symsize=4, color='blue'

    djs_oplot, weff, 1D17*ugriz_flam, psym=6, symsize=4, color='yellow'
    djs_oplot, weff, 1D17*keck_ugriz_flam, psym=7, symsize=4, color='dark red'
    djs_oplot, weff, 1D17*newkeck_ugriz_flam, psym=7, symsize=4, color='orange'

;   djs_oplot, keck.wavelength, 1D17*poly(keck.wavelength,sdsscoeff), line=0, color='purple'
;   djs_oplot, keck.wavelength, 1D17*poly(keck.wavelength,keckcoeff), line=0, color='red'
;   djs_oplot, keck.wavelength, medsmooth(keck.flux*1D17,10)*keckcorr, psym=10, color='purple'
    
;   djs_oplot, keck.wavelength, 1D17*keckc, line=0, color='purple'
;   djs_oplot, mmt_a.wave, 1D17*mmtc, line=0, color='red'
    
    info = im_filterspecs(filterlist=sdss_filterlist())
    for ii = 0, 4 do begin
       ff = info[ii].filtf[0:info[ii].filtn-1]
       djs_oplot, info[ii].filtw[0:info[ii].filtn-1], ff/max(ff)*!y.crange[1]*0.2
    endfor
    im_plotconfig, psfile=psfile, /psclose, /pdf
    
return
end
    
