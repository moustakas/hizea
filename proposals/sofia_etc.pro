pro sofia_etc
; jm15may19siena - estimate the exposure time to observe the highzea sample with
; the Sofia/HAWC instrument

    light = 2.99792458D18       ; speed of light [A/s]

    masspath = massprofiles_path()
    sofiapath = getenv('HIZEA_DATA')+'/sofia/'
    info = mrdfits(masspath+'massprofiles_lir.fits.gz',1)
    phot = mrdfits(masspath+'massprofiles_photometry.fits.gz',1)
    seds = mrdfits(masspath+'massprofiles_lir_seds.fits.gz',1)

    keep = where(phot.maggies[11] gt 0)
    info = info[keep]
    phot = phot[keep]
    seds = seds[keep]
    ngal = n_elements(info)

    pred = replicate({galaxy: '', z: 0.0, hawc_chary: fltarr(5), hawc_dale: fltarr(5)},ngal)
    pred.galaxy = info.galaxy
    pred.z = info.z

    bands = ['A','B','C','D','E']
    filters = 'sofia_hawc_band'+bands+'.par'
    weff_hawc = k_lambda_eff(filterlist=filters)

    weff = k_lambda_eff(filterlist=[hizea_filterlist(),wise_filterlist()])
    normband = 12 ; 22-micron band

; read the QSO SED
    light = im_light(/ang)
    readcol, getenv('HIZEA_DIR')+'/etc/qso2_sed.txt', qlam, nuLnu, $
      skip=4, format='F,D', /silent
    qflux = nuLnu*qlam/im_light(/micron) ; [nu*L_nu --> erg/s/Hz]

    psfile = sofiapath+'sofia_etc.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.2, $
      xmargin=[1.6,0.3], width=6.6, thick=8
    for ii = 0, ngal-1 do begin
       wave_chary = seds[ii].modelwave_chary
       fnu_chary = 10^(-0.4*(seds[ii].modelmab_chary+48.6))
       maggies_chary = reform(k_project_filters(k_lambda_to_edges(wave_chary),$
         fnu_chary*light/wave_chary^2,filterlist=filters))*10^(-0.4*48.6)*1D26 ; [mJy]

       wave_dale = seds[ii].modelwave_dale
       fnu_dale = 10^(-0.4*(seds[ii].modelmab_dale+48.6))
       maggies_dale = reform(k_project_filters(k_lambda_to_edges(wave_dale),$
         fnu_dale*light/wave_dale^2,filterlist=filters))*10^(-0.4*48.6)*1D26 ; [mJy]

       pred[ii].hawc_chary = maggies_chary
       pred[ii].hawc_dale = maggies_dale

       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         /xlog, /ylog, xrange=[1.0,5E3], yrange=[0.01,500], $
         xtitle='Wavelength (\mu'+'m)', ytitle='F_{\nu} (mJy)'
       djs_oplot, wave_chary/1D4, fnu_chary*1D26
       djs_oplot, wave_dale/1D4, fnu_dale*1D26, line=5

       djs_oplot, weff/1D4, phot[ii].maggies*10^(-0.4*48.6)*1D26, psym=6, symsize=3, color='blue'
       djs_oplot, weff_hawc/1D4, maggies_chary, psym=6, symsize=3, color='orange'
       djs_oplot, weff_hawc/1D4, maggies_dale, psym=cgsymcat(9), $
         symsize=3, color=cgcolor('forest green')

       zqlam = (1+pred[ii].z)*qlam
       qnorm = phot[ii].maggies[normband]*10^(-0.4*48.6)*1D26/interpol(qflux,zqlam,weff[normband]/1D4)
       djs_oplot, zqlam, qflux*qnorm, line=3, color=cgcolor('firebrick')
       
       im_legend, info[ii].galaxy, /left, /top, box=0
       im_legend, ['Chary & Elbaz','Dale & Helou'], /right, /bottom, $
         box=0, line=[0,5], pspacing=1.7, charsize=1.5

;      cc = get_kbrd(1)
    endfor
    im_plotconfig, psfile=psfile, /psclose, /pdf

    mwrfits, pred, sofiapath+'sofia_etc.fits', /create
    struct_print, pred


stop

return
end
