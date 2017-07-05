pro hawc_cycle5
; jm15may19siena - estimate the exposure time to observe the highzea sample with
; the Sofia/HAWC instrument

    normband = 12 ; 22-micron band
    wiseband = [11,12]

    light = 2.99792458D18       ; speed of light [A/s]
    maggies2mJy = 10^(-0.4*48.6)*1D26

    masspath = massprofiles_path()
    hawcpath = getenv('HIZEA_DATA')+'/hawc_cycle5/'
    phot = mrdfits(masspath+'massprofiles_photometry.fits.gz',1)

    match, phot.galaxy, ['J1613+2834','J1506+5402'], m1, m2
;   match, phot.galaxy, ['J1613+2834','J1506+5402','J1219+0336','J1341-0321','J0905+5759'], m1, m2
    phot = phot[m1]
    ngal = n_elements(phot)

    srt = reverse(sort(phot.maggies[normband]))
    phot = phot[srt]

; read the IR SEDs
    readcol, getenv('HIZEA_DIR')+'/etc/ir_seds.txt', irwave, fnu_ms, $
      fnu_sb, fnu_agn, format='F,F,F', skipline=5, /silent
    irwave *= 1D4 ; [Angstrom]

    pred = replicate({galaxy: '', z: 0.0, hawc_ms: fltarr(5), $
      hawc_sb: fltarr(5), hawc_agn: fltarr(5), $
      hawc_ms_snr: fltarr(5), hawc_sb_snr: fltarr(5), $
      hawc_agn_snr: fltarr(5)},ngal)
    pred.galaxy = phot.galaxy
    pred.z = phot.z

;   bands = ['C','E']
    bands = ['A','B','C','D','E']
    filters = 'sofia_hawc_band'+bands+'.par'
    weff_hawc = k_lambda_eff(filterlist=filters)

    ace = [0,2,4] ; filter A,C,E
    ce = [2,4] ; filter C,E

    weff = k_lambda_eff(filterlist=[hizea_filterlist(),wise_filterlist()])

    psfile = hawcpath+'hawc_cycle5.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=4.5, $
      xmargin=[1.6,0.3], width=6.6, thick=8
    for ii = 0, ngal-1 do begin

; put Christy's S/N estimates in (just for bands A, C, E)
       case pred[ii].galaxy of
          'J1613+2834': begin
             pred[ii].hawc_ms_snr = [3.8,0.0,12.0,0.0,38.2]
             pred[ii].hawc_sb_snr = [11.8,0.0,  28.9,0.0,  34.6]
             pred[ii].hawc_agn_snr = [4.0,0.0,4.0,0.0,4.0]
          end
          'J1506+5402': begin
             pred[ii].hawc_ms_snr = [3.1,0.0,   8.5,0.0,  32.4]
             pred[ii].hawc_sb_snr = [10.6,0.0,  28.9,0.0,  40.0]
             pred[ii].hawc_agn_snr = [4.0,0.0,4.0,0.0,4.0]
          end
          'J1219+0336': begin
             pred[ii].hawc_ms_snr = [3.8,0.0,  12.0,0.0,  38.3]
             pred[ii].hawc_sb_snr = [11.8,0.0,  29.0,0.0,  34.8]
             pred[ii].hawc_agn_snr = [4.0,0.0,4.0,0.0,4.0]
          end
          'J1341-0321': begin
             pred[ii].hawc_ms_snr = [2.6,0.0,   6.8,0.0,  26.9]
             pred[ii].hawc_sb_snr = [9.1,0.0,  25.5,0.0,  36.9]
             pred[ii].hawc_agn_snr = [4.0,0.0,4.0,0.0,4.0]
          end
          'J0905+5759': begin
             pred[ii].hawc_ms_snr = [1.8,0.0,   4.7,0.0,  19.5]
             pred[ii].hawc_sb_snr = [6.9,0.0,  20.2,0.0,  30.6]
             pred[ii].hawc_agn_snr = [4.0,0.0,4.0,0.0,4.0]
          end
          else: 
       endcase

; normalize at 22 microns
       zirwave = (1+pred[ii].z)*irwave
       zirwave_edges = k_lambda_to_edges(zirwave)

       fnu2flam = 1D-3*1D-23*light/zirwave^2 ; [mJy --> erg/s/cm2/A]
       flam_ms = fnu_ms*fnu2flam
       flam_sb = fnu_sb*fnu2flam
       flam_agn = fnu_agn*fnu2flam

       norm_ms = phot[ii].maggies[normband]*maggies2mJy/$
         interpol(fnu_ms,zirwave,weff[normband])
       norm_sb = phot[ii].maggies[normband]*maggies2mJy/$
         interpol(fnu_sb,zirwave,weff[normband])
       norm_agn = phot[ii].maggies[normband]*maggies2mJy/$
         interpol(fnu_agn,zirwave,weff[normband])

       maggies_ms = reform(k_project_filters(zirwave_edges,$
         norm_ms*flam_ms,filterlist=filters))*maggies2mJy   ; [mJy]
       maggies_sb = reform(k_project_filters(zirwave_edges,$
         norm_sb*flam_sb,filterlist=filters))*maggies2mJy   ; [mJy]
       maggies_agn = reform(k_project_filters(zirwave_edges,$
         norm_agn*flam_agn,filterlist=filters))*maggies2mJy ; [mJy]

       pred[ii].hawc_ms = maggies_ms
       pred[ii].hawc_sb = maggies_sb
       pred[ii].hawc_agn = maggies_agn

       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         /xlog, /ylog, xrange=[3.0,3E3], yrange=[0.1,900], $
         xtitle='Wavelength (\mu'+'m)', ytitle='F_{\nu} (mJy)'
       djs_oplot, zirwave/1D4, norm_ms*fnu_ms
       djs_oplot, zirwave/1D4, norm_sb*fnu_sb, line=3
       djs_oplot, zirwave/1D4, norm_agn*fnu_agn, line=5

       djs_oplot, weff/1D4, phot[ii].maggies*maggies2mJy, psym=6, $
         symsize=3, color=cgcolor('orange')
       djs_oplot, weff_hawc/1D4, maggies_ms, psym=symcat(6), $
         symsize=3, color=cgcolor('forest green')
       djs_oplot, weff_hawc/1D4, maggies_sb, psym=symcat(5), $
         symsize=3, color=cgcolor('dodger blue')
       djs_oplot, weff_hawc/1D4, maggies_agn, psym=symcat(9), $
         symsize=3, color=cgcolor('firebrick')
;      djs_oplot, weff_hawc/1D4, maggies_dale, psym=cgsymcat(9), $
;        symsize=3, color=cgcolor('forest green')

       im_legend, phot[ii].galaxy, /left, /top, box=0
       im_legend, ['Main seq. SF','Starburst','AGN'], /right, /top, $
         box=0, line=[0,3,5], pspacing=1.9, charsize=1.5, $
         position=[0.75,0.45], /norm, thick=8, $
         color=['forest green','dodger blue','firebrick']

    endfor
    im_plotconfig, psfile=psfile, /psclose, /pdf

    mwrfits, pred, hawcpath+'hawc_cycle5.fits', /create
;   struct_print, pred

; make a figure for the proposal using Christy's S/N estimates
    ynorm = 1 ; [mJy]
;   ynorm = 1E3 ; [mJy --> Jy]
    symsize1 = 1.4

    psfile = hawcpath+'hawc_cycle5_etc.ps'
    im_plotconfig, 7, pos, psfile=psfile, charsize=1.5, $
      height=2.3, xmargin=[1.2,0.3]
    for ii = 0, ngal-1 do begin

; normalize at 22 microns
       zirwave = (1+pred[ii].z)*irwave
       zirwave_edges = k_lambda_to_edges(zirwave)

       fnu2flam = 1D-3*1D-23*light/zirwave^2 ; [mJy --> erg/s/cm2/A]
       flam_ms = fnu_ms*fnu2flam
       flam_sb = fnu_sb*fnu2flam
       flam_agn = fnu_agn*fnu2flam

       norm_ms = phot[ii].maggies[normband]*maggies2mJy/$
         interpol(fnu_ms,zirwave,weff[normband])
       norm_sb = phot[ii].maggies[normband]*maggies2mJy/$
         interpol(fnu_sb,zirwave,weff[normband])
       norm_agn = phot[ii].maggies[normband]*maggies2mJy/$
         interpol(fnu_agn,zirwave,weff[normband])

; make the plot       
;      if ii le 2 then xtickname = replicate(' ',10) else delvarx, xtickname
       if odd(ii) then ytickname = replicate(' ',10) else delvarx, ytickname
       if ii ge 0 then xtitle = 'Wavelength (\mu'+'m)' else xtitle = ''
       if ii eq 0 then ytitle = 'Flux Density (mJy)' else ytitle = ''
;      if ii eq 2 then ytitle = 'F_{\nu} (mJy)' else ytitle = ''

       djs_plot, [0], [0], /nodata, position=pos[*,ii], xsty=1, ysty=1, $
         /xlog, /ylog, noerase=ii gt 0, xrange=[7.0,1E3], yrange=[0.15,1000]/ynorm, $
         xtickname=xtickname, ytickname=ytickname, xtitle=xtitle, ytitle=ytitle
       djs_oplot, zirwave/1D4, norm_ms*fnu_ms/ynorm, color=cgcolor('olive'), thick=7
       djs_oplot, zirwave/1D4, norm_sb*fnu_sb/ynorm, line=3, color=cgcolor('dodger blue'), thick=7
       djs_oplot, zirwave/1D4, norm_agn*fnu_agn/ynorm, line=5, color=cgcolor('tomato'), thick=7

       djs_oplot, weff/1D4, phot[ii].maggies*maggies2mJy/ynorm, psym=symcat(6,thick=6), $
         symsize=1.5

       these = ce 
;      if strmatch(phot[ii].galaxy,'J0905*') then these = ce else these = ace

       oploterror, weff_hawc[these]/1D4, pred[ii].hawc_ms[these]/ynorm, $
         pred[ii].hawc_ms[these]/pred[ii].hawc_ms_snr[these]/ynorm, $
         psym=symcat(16), symsize=symsize1, color=cgcolor('forest green'), $
         errcolor=cgcolor('forest green')

       oploterror, weff_hawc[these]/1D4, pred[ii].hawc_sb[these]/ynorm, $
         pred[ii].hawc_sb[these]/pred[ii].hawc_sb_snr[these]/ynorm, $
         psym=symcat(14), symsize=symsize1*1.2, color=cgcolor('navy'), $
         errcolor=cgcolor('navy')

       oploterror, weff_hawc[these]/1D4, pred[ii].hawc_agn[these]/ynorm, $
         pred[ii].hawc_agn[these]/pred[ii].hawc_agn_snr[these]/ynorm, $
         psym=symcat(15), symsize=symsize1, color=cgcolor('firebrick'), $
         errcolor=cgcolor('firebrick')

       im_legend, phot[ii].galaxy, /left, /top, box=0, position=[pos[0,ii]-0.005,pos[3,ii]-0.01], $
         charsize=1.3, /norm

; add some tick marks for WISE and HAWC
       if ii eq 0 then begin
;         xyouts, 15, 0.5, 'WISE!cW1,W2', /data, charsize=1.3, align=0.5
          xyouts, 17, 0.3, 'WISE', /data, charsize=1.3, align=0.5
          xyouts, 12, 0.8, 'W3', /data, charsize=1.3, align=0.5
          xyouts, 22, 3.0, 'W4', /data, charsize=1.3, align=0.5

;         xyouts, 100, 12.0, 'HAWC!cBands ACE', /data, charsize=1.3, align=0.5
          xyouts, 100, 6.0, 'HAWC', /data, charsize=1.3, align=0.5
;         xyouts, 54, 17.0, 'A', /data, charsize=1.3, align=0.5
          xyouts, 88, 19.0, 'C', /data, charsize=1.3, align=0.5
          xyouts, 215, 12.0, 'E', /data, charsize=1.3, align=0.5
       endif
    endfor    

    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,ii], $
      xsty=5, ysty=5
    im_legend, ['Main seq. SF','Starburst','AGN'], $
      box=0, line=[0,3,5], pspacing=2.1, charsize=1.0, $
      position=[0.68,0.52], /norm, thick=8, $
      color=cgcolor(['olive','dodger blue','tomato'])

    xx = 0.725 & yy = 0.49
    plots, xx, yy, psym=symcat(16), /norm, $
      color=cgcolor('forest green'), symsize=1.3
    plots, xx, yy-0.045, psym=symcat(14), /norm, $
      color=cgcolor('navy'), symsize=1.3*1.3
    plots, xx, yy-0.095, psym=symcat(15), /norm, $
      color=cgcolor('firebrick'), symsize=1.3

    im_plotconfig, psfile=psfile, /psclose, /pdf




stop

return
end
