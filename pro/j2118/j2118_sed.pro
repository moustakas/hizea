pro j2118_sed, final=final
; jm19feb07siena - estimate the 2mm flux for J2118

    light = 2.99792458D18       ; speed of light [A/s]
    maggies2mJy = 10^(-0.4*48.6)*1D26

    masspath = massprofiles_path()
    outpath = getenv('HIZEA_DATA')+'/j2118_sed/'
    phot = mrdfits(masspath+'massprofiles_photometry_all.fits.gz',1)

;   match, phot.galaxy, ['J1613+2834','J1506+5402'], m1, m2
;   match, phot.galaxy, ['J1613+2834','J1506+5402','J1219+0336','J1341-0321','J0905+5759'], m1, m2

    match, phot.galaxy, ['J2118+0017'], m1, m2
    phot = phot[m1]
    ngal = n_elements(phot)

;   normband = 12 ; 22-micron band
;   srt = reverse(sort(phot.maggies[normband]))
;   phot = phot[srt]

; read the IR SEDs
    readcol, getenv('HIZEA_DIR')+'/etc/ir_seds.txt', irwave, fnu_ms, $
      fnu_sb, fnu_agn, format='F,F,F', skipline=5, /silent
    irwave *= 1D4 ; [Angstrom]

    pred = replicate({galaxy: '', z: 0.0, ms: 0.0, $
      sb: 0.0, agn: 0.0},ngal)
    pred.galaxy = phot.galaxy
    pred.z = phot.z

    weff_2mm = 2E-3/1E-10 ; [Angstrom]
    weff = k_lambda_eff(filterlist=[hizea_filterlist(),wise_filterlist()])

    psfile = outpath+'j2118_sed.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=4.3, $
      xmargin=[1.6,0.3], width=6.6, thick=8
    for ii = 0, ngal-1 do begin

; normalize at 22 microns
       zirwave = (1+pred[ii].z)*irwave
       zirwave_edges = k_lambda_to_edges(zirwave)

       fnu2flam = 1D-3*1D-23*light/zirwave^2 ; [mJy --> erg/s/cm2/A]
       flam_ms = fnu_ms*fnu2flam
       flam_sb = fnu_sb*fnu2flam
       flam_agn = fnu_agn*fnu2flam

       if phot[ii].maggies[12] gt 0.0 then normband = 12 else normband = 11 ; 22-micron band, else 12 micron
       
       norm_ms = phot[ii].maggies[normband]*maggies2mJy/$
         interpol(fnu_ms,zirwave,weff[normband])
       norm_sb = phot[ii].maggies[normband]*maggies2mJy/$
         interpol(fnu_sb,zirwave,weff[normband])
       norm_agn = phot[ii].maggies[normband]*maggies2mJy/$
         interpol(fnu_agn,zirwave,weff[normband])

       pred[ii].ms = interpol(fnu_ms*norm_ms,zirwave,weff_2mm)
       pred[ii].sb = interpol(fnu_sb*norm_sb,zirwave,weff_2mm)
       pred[ii].agn = interpol(fnu_agn*norm_agn,zirwave,weff_2mm)

       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         /xlog, /ylog, xrange=[3.0,5E3], yrange=[0.01,900], $
         xtitle='Observed-frame Wavelength (\mu'+'m)', ytitle='F_{\nu} (mJy)'; $
;        xticklen=1.0, yticklen=1.0, xgridstyle=1, ygridstyle=1, $

       djs_oplot, zirwave/1D4, norm_ms*fnu_ms
       djs_oplot, zirwave/1D4, norm_sb*fnu_sb, line=3
       if not keyword_set(final) then $
         djs_oplot, zirwave/1D4, norm_agn*fnu_agn, line=5

       djs_oplot, weff/1D4, phot[ii].maggies*maggies2mJy, psym=6, $
         symsize=3, color=cgcolor('orange')
       djs_oplot, [weff_2mm/1D4], [pred[ii].ms], psym=symcat(9), $
         symsize=3, color=cgcolor('forest green')
       djs_oplot, [weff_2mm/1D4], [pred[ii].sb], psym=symcat(5), $
         symsize=3, color=cgcolor('dodger blue')
       if not keyword_set(final) then begin
          djs_oplot, [weff_2mm/1D4], [pred[ii].agn], psym=symcat(9), $
            symsize=3, color=cgcolor('firebrick')
       endif

       if keyword_set(final) then begin
          xyouts, 800, 200, 'Band 9!cPrediction', /data, align=0.5, $
            charsize=1.3
          im_legend, 'WISE W1-W4', color='orange', /left, /top, box=0, $
            psym=6, symsize=3, charsize=1.7
          im_legend, ['Star-forming','Starburst'], /right, /top, $
            box=0, line=[0,3], pspacing=1.9, charsize=1.5, $
            position=[0.8,0.45], /norm, thick=8;, $
;           color=['forest green','dodger blue']
       endif else begin
          im_legend, phot[ii].galaxy, /left, /top, box=0
          im_legend, ['Main seq. SF','Starburst','AGN'], /right, /top, $
            box=0, line=[0,3,5], pspacing=1.9, charsize=1.5, $
            position=[0.75,0.45], /norm, thick=8, $
            color=['forest green','dodger blue','firebrick']
       endelse

    endfor
    im_plotconfig, psfile=psfile, /psclose, /pdf

    mwrfits, pred, outpath+'j2118_sed.fits', /create

stop

return
end
