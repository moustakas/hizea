pro sigmasfr_datafigure
; ad12mar19ucsd - written to make Figure 1 for letter

  !p.thick=3
  !x.thick=3
  !y.thick=3
  !p.charthick=3

  path = getenv('HIZEA_DATA')+'/sigmasfr/'
  phot = mrdfits(path+'sigmasfr_photometry.fits.gz', 1)
  hst = rsex(path+'hst_sample.dat')
  kcorr = mrdfits(path+'sigmasfr_kcorrect.fits.gz',1)
  lir = mrdfits(path+'sigmasfr_lir.fits.gz',1)

  filters = sigmasfr_filterlist()
  weff = k_lambda_eff(filterlist=filters)

    red, omega0=0.3, omegalambda=0.7, h100=0.7;, /verb


  asize = 1.3

  dfpsplot, path+'sigmasfr_datafigure.ps', /color, /landscape

  erase

  for jj=0L, 2L do begin;n_elements(phot)-1L do begin

    used = where((phot[jj].maggies gt 0.0) and $ 
      (phot[jj].ivarmaggies gt 0.0),nused)

    mab =  maggies2mag(phot[jj].maggies[used],$
      ivar=phot[jj].ivarmaggies[used],magerr=mab_err)

    xtitle1 = textoidl('Observed-Frame Wavelength (\AA)')
    xtitle2 = textoidl('Rest-Frame Wavelength (\mum)')
    ytitle1 = textoidl('Apparent AB Magnitude')
    ytitle2 = textoidl('F_{\nu} [mJy]')
    ytitle3 = textoidl('F_{\lambda}')

    yrange = [24.0,13.0]
    yrange2 = 10.^(yrange*(-1.)/2.5)*3631.d3
    fac = 1.
    ;fac = 10000.

    ticks1 = [1000,2000,4000,8000,16000,30000,60000,120000,240000]/fac
    ticks2 = [1000,2000,4000,8000,16000,30000,60000,120000,240000]/fac

    ;if jj eq 0L then pos=[0.1, 0.5, 0.5, 0.9]
    pos = subplot(2,2,[jj], egap=0.07)

    xrange1 = [1200,300000]/fac
    ;xrange2 = xrange1/(1.0+phot[jj].z)
    xrange2 = [700., 300000.]
    plot, [0], [0], /nodata, xrange=xrange2, yrange=yrange, $
      xsty=5, ysty=5, xlog=1, $ 
      position=pos, xticks=n_elements(ticks1)-1, $
      xtickv=ticks1, xtickname=replicate(' ',8), $
      ytickname=replicate(' ',8), $
      color=djs_icolor('black'), charsize=asize, $
      xmargin=[7.,2.], ymargin=[4.,3.], /noerase
; take care of x-axes
    if jj eq 0L or jj eq 1L then begin
    axis, xaxis=1, xsty=1, xtitle=xtitle2, xrange=xrange2, xlog=1, $
      xticks=n_elements(ticks2)-1, xtickv=ticks2, $;xtickformat='(I0)', $
      charsize=asize, xtickname=['0.1', '0.2', '0.4', '0.8', '1.6', '3', $
        '6', '12', '24']
    axis, xaxis=0, xsty=1, xrange=xrange2, xlog=1, $
      xticks=n_elements(ticks1)-1, $
      xtickv=ticks1, xtickname=replicate(' ',9)
    endif else begin
    axis, xaxis=0, xsty=1, xtitle=xtitle2, xrange=xrange2, xlog=1, $
      xticks=n_elements(ticks2)-1, xtickv=ticks2, $;xtickformat='(I0)', $
      charsize=asize, xtickname=['0.1', '0.2', '0.4', '0.8', '1.6', '3', $
        '6', '12', '24']
    axis, xaxis=1, xsty=1, xrange=xrange2, xlog=1, $
      xticks=n_elements(ticks1)-1, $
      xtickv=ticks1, xtickname=replicate(' ',9)
    endelse
; take care of y-axes
    if jj eq 0L or jj eq 2L then begin
      axis, yaxis=0, ysty=1, ytitle=ytitle1, yrange=yrange, charsize=asize
      axis, yaxis=1, ysty=1, yrange=yrange, ytickname=replicate(' ',9)
    endif else begin
      axis, yaxis=0, ysty=1, yrange=yrange, ytickname=replicate(' ',9)
      axis, yaxis=1, ysty=1, ytitle=ytitle2, yrange=yrange2, /ylog, $
      ytickname=['0.001', '0.01', '0.1', '1', '10'], charsize=asize
    endelse

    oploterror, weff[used]/fac/(1.+phot[jj].z), mab, mab_err, psym=symcat(16), $
      symsize=1.5, color=djs_icolor('red'), $
      errcolor=djs_icolor('red'), errthick=8     
;stop
   ; xyouts, 840., 17.0, 'HST WFC3/F814W', size=0.8
   ; ;xyouts, mean([posb[0],posb[2]]), posb[1]-dy, $
   ; xyouts, 4200., 17.0, 'GALFIT model', size=0.8
   ; ;xyouts, mean([posc[0],posc[2]]), posc[1]-dy, $
   ; xyouts, 20000., 17.0, 'residual', size=0.8
   ; ;  size=1.5, align=0.5, /normal

;if max(weff[used]) gt 20000. then begin

; plot M82 SED scaled to WISE ch2
  readcol, path+'/swire/M82_template_norm.sed', m82wave, m82flambda
  m82mag = -2.5*alog10(m82flambda * (m82wave)^2)
  sort = sort(abs(m82wave - weff[8]/(1.+phot[jj].z)))
  loc = where(weff[used] eq weff[8])
  offset = mab[loc] - m82mag[sort[0]]
  good = where(m82wave gt 25000., mgood)
  oplot, m82wave[good], m82mag[good]+offset[0], color=djs_icolor('blue'), thick=2

; plot ARP220 SED scaled to WISE ch2
  readcol, path+'/swire/ARP220_template_norm.sed', arp220wave, arp220flambda
  arp220mag = -2.5*alog10(arp220flambda * (arp220wave)^2)
  sort = sort(abs(arp220wave - weff[8]/(1.+phot[jj].z)))
  loc = where(weff[used] eq weff[8])
  offset = mab[loc] - arp220mag[sort[0]]
  good = where(arp220wave gt 25000., mgood)
  oplot, arp220wave[good], arp220mag[good]+offset[0], color=djs_icolor('dark green'), linestyle=2, thick=6

; plot QSO2 SED scaled to WISE ch2
  readcol, path+'/swire/torus_template_norm.sed', toruswave, torusflambda
  torusmag = -2.5*alog10(torusflambda * (toruswave)^2)
  sort = sort(abs(toruswave - weff[8]/(1.+phot[jj].z)))
  loc = where(weff[used] eq weff[8])
  offset = mab[loc] - torusmag[sort[0]]
  good = where(toruswave gt 25000., mgood)
  oplot, toruswave[good], torusmag[good]+offset[0], color=djs_icolor('red'), $
    linestyle=1, thick=6

;stop
;endif

if phot[jj].galaxy eq 'J1506+5402' then $
  legend, ['M82 starburst', 'Arp 220 starburst', 'obscured quasar'], $
    linestyle=[0, 2, 1], color=[djs_icolor('blue'), djs_icolor('dark green'), djs_icolor('red')], $
    charsize=1.3, box=0, pos=[3500., 19.5], thick=6

; put galaxy name on plot
    ;xyouts, 950., 18.0, phot[jj].galaxy, charsize=1.2
    ;xyouts, 13000., 15.0, phot[jj].galaxy, charsize=1.5
    ;xyouts, 5000., 20.0, phot[jj].galaxy+' (z='+string(phot[jj].z, format='(f5.3)')+')', charsize=1.3
    ;xyouts, 5000., 21.0, textoidl('\Sigma_{SFR}=')+string(round(hst[jj].wise_sfr / (2.*!pi*hst[jj].r_e^2)/ 10.)*10., $
    ;  format='(i4)')+' M'+sunsymbol()+textoidl(' yr^{-1} kpc^{-2}'), charsize=1.3
    ;xyouts, 5000., 22.5, textoidl('v_{out}=')+string(round(hst[jj].vout/10.)*10.,format='(i5)')+$
    ;  textoidl(' km s^{-1}'), charsize=1.5
    xyouts, 1300., 23.3, phot[jj].galaxy+', '+textoidl('\Sigma_{SFR}=')+$
      string(round(hst[jj].wise_sfr / (2.*!pi*hst[jj].r_e^2)/ 10.)*10., format='(i4)')+$
      ' M'+sunsymbol()+textoidl(' yr^{-1} kpc^{-2}'), charsize=1.3


; put galfit results on plot
  ; horizonal spacing on the top left
    posa = [0.03, 0.68, 0.26, 0.96]
    posb = [0.26, 0.68, 0.49, 0.96] 
    posc = [0.49, 0.68, 0.72, 0.96]

    posone = [0.05, 0.55, 0.45, 0.96]
    ;posone = [0.3, 0.05, 0.7, 0.46]

    if phot[jj].galaxy eq 'J0905+5759' then $
       plot_sigmasfr_galfit, path+'j0905+5759_drz_sci.fits', 5955., 5880., $
          npix=120., posa=subsubplot(pos, posone), /onlydata, /ext0
    if phot[jj].galaxy eq 'J1341-0321' then $
      ;plot_sigmasfr_galfit, path+'imgblock_j1341-0321_sersic_neq4.fits', $
      ;  105., 106., npix=60, posa=subsubplot(pos, posone), /onlydata 
       plot_sigmasfr_galfit, path+'j1341-0321_drz_sci.fits', 5258., 5457., $
          npix=120., posa=subsubplot(pos, posone), /onlydata, /ext0
    if phot[jj].galaxy eq 'J1506+5402' then $
      ;plot_sigmasfr_galfit, path+'imgblock_j1506+5402_sersic_neq4.fits', $
      ;  99., 103., npix=60., posa=subsubplot(pos, posone), /onlydata 
       plot_sigmasfr_galfit, path+'j1506+5402_drz_sci.fits', 5684., 5524., $
          npix=120., posa=subsubplot(pos, posone), /onlydata, /ext0 
    if phot[jj].galaxy eq 'J0944+0930' then $
      plot_sigmasfr_galfit, path+'imgblock_j0944+0930_sersic_neq4.fits', $
        99., 103., posa=subsubplot(pos, posone), /onlydata 

kpcarc = dangular(0.7116, /kpc)/206265.d

    hst_xrange = [-120.5*kpcarc*0.02, 120.5*kpcarc*0.02]
    hst_yrange = [-120.5*kpcarc*0.02, 120.5*kpcarc*0.02]
    plot, [0], [0], pos=subsubplot(pos,posone), /noerase, $
      color=djs_icolor('black'), xsty=5, ysty=5, $
      xrange=hst_xrange, $
      yrange=hst_yrange
    plots, !x.crange, [hst_yrange[0], hst_yrange[0]]
    plots, !x.crange, [hst_yrange[1], hst_yrange[1]]
    plots, [hst_xrange[0], hst_xrange[0]], !y.crange
    plots, [hst_xrange[1], hst_xrange[1]], !y.crange
    plots, [5., 15.], $
      [hst_yrange[0]*0.75, hst_yrange[0]*0.75], $
      color=djs_icolor('black'), thick=10
    xyouts, 10., hst_yrange[0]*0.65, '10 kpc', charsize=1.1, $
      align=0.5, charthick=5.

    ;  xrange=[-120.5*kpcarc*0.02, 120.5*kpcarc*0.02], $
    ;  yrange=[-120.5*0.02, 120.5*0.02]
      ;xtickname=replicate(' ',6), xticks=1, xminor=0, $
      ;ytickname=replicate(' ',6), yticks=1., yminor=0

    ;plot, [0], [0], pos=subsubplot(pos,posb), /noerase, $
    ;  color=djs_icolor('black'), $
    ;  xtickname=replicate(' ',6), $
    ;  ytickname=replicate(' ',6)
    ;
    ;plot, [0], [0], pos=subsubplot(pos,posc), /noerase, $
    ;  color=djs_icolor('black'), $
    ;  xtickname=replicate(' ',6), $
    ;  ytickname=replicate(' ',6)


;; put lores spectrum on plot
;  ; spectra across the bottom
;    pos1 = [0.20, 0.04, 0.61, 0.33]
;    pos2 = [0.61, 0.04, 0.78, 0.33] 
;    pos3 = [0.78, 0.04, 0.95, 0.33]
;
;    posspec = [0.20, 0.04, 0.95, 0.33]
;
;    ;xrange1 = [4000,8700]
;    ;xrange2 = xrange1/(1.0+phot[jj].z)
;    xrange2 = [2600., 5200.]
;    yrange = [0.0,20]
;    xtitle1 = textoidl('Observed-Frame Wavelength (\AA)')
;    xtitle2 = textoidl('Rest-Frame Wavelength (\AA)')
;    xtitle3 = textoidl('Rest \lambda [')+angstrom()+']'
;    ytitle1 = textoidl('F_{\lambda}')
;    xtlen=0.05
;    ym1 = 5.
;    ytint1 = 5
;    ytlen1 = 0.04
;    xtint2 = 1000.
;    xm2 = 10.
;
;    plot, [0], [0], /nodata, xrange=xrange2, yrange=yrange, $
;      xsty=5, ysty=1, xtitle=xtitle1, ytitle=ytitle1, $
;      yminor=ym1, ytickint=ytint1, yticklen=ytlen1, $ 
;      position=subsubplot(pos,posspec), $
;      /noerase, color=djs_icolor('black'), xticklen=xtlen, charsize=asize
;    axis, /xaxis, xsty=1, xtitle=xtitle3, xrange=xrange2, xticklen=xtlen, $
;      xtickint=xtint2, xminor=xm2, charsize=0.8
;    axis, xaxis=0, xsty=1, xrange=xrange2, xtickname=replicate(' ',6), $
;        xticklen=xtlen    
;
;    if phot[jj].galaxy eq 'J0905+5759' then begin
;      j0905path = j0905_path()
;      lris = mrdfits(j0905path+'spectra/j0905+5759_lris_flux_v120216.fits',1)
;      scale = 1D17
;      djs_oplot, lris.wavelength/(1.+phot[jj].z), smooth(scale*lris.flux,1), color='grey'       
;    endif
;
;    if phot[jj].galaxy eq 'J1506+5402' then begin
;      mmt = rsex(path+'j1506+5402.txt')
;      djs_oplot, mmt.wave/(1.+phot[jj].z), smooth(mmt.flux,1), color='grey'
;    endif
;
;    if phot[jj].galaxy eq 'J1341-0321' then begin
;      mage = rsex(path+'j1341-0321.txt')
;      djs_oplot, mage.wave/(1.+phot[jj].z), smooth(mage.flux,3), color='grey'
;    endif
;
;  ; plot hires spectrum    
;    xrange3 = [2765., 2785.]
;    yrange3 = [0., 1.2]
;    xminor3 = 4.
;    xtint3 = 10.
;    xtitle4 = textoidl('\lambda [')+angstrom()+']'
;    xtlen=0.04
;    yminor3 = 5
;
;    plot, [0], [0], /nodata, pos=subsubplot(pos,pos2), /noerase, $
;        xrange=xrange3, xsty=1, xminor=xminor3, xtickint=xtint3, $
;        xtickname=replicate(' ',7), xticklen=xtlen, $
;        ytickname=replicate(' ',7), $
;        yrange=yrange3, ysty=5, yminor=yminor3, ytickint=ytint
;      axis, yaxis=1, yrange=yrange3, ysty=1, yminor=yminor3, ytickint=ytint, $
;        ytickname=replicate(' ',7), yticklen=ytlen
;      axis, xaxis=1, xrange=xrange3, xsty=1, xminor=xminor3, $
;        xtickint=xtint3, xticklen=xtlen, xtitle=xtitle4, charsize=0.8
;      
;;xticks=n_elements(ticks2)-1, xtickv=ticks2, $;xtickformat='(I0)', $
;;      charsize=asize, xtickname=['0.1', '0.2', '0.4', '0.8', '1.6', '3', $
;;        '6', '12', '24']
;
;    if phot[jj].galaxy eq 'J0905+5759' then begin
;      z_bg = 0.7116d    
;      hires = rsex(path+'J0905a.dat')
;      oplot, hires.wave/(1.+z_bg), medsmooth(hires.flux, 5), psym=10
;    endif 
;
;  ; plot hires spectrum in velocity space
;    xrange4 = [-3000., -2200.]
;    xminor4 = 4.
;    xtint4 = 400.
;    ;xtitle4 = textoidl('Mg II Velocity [km s^{-1}]')
;    xtitle5 = textoidl('v [km s^{-1}]')
;    ytint = 0.5
;    ytm = 5
;    xtlen = 0.04
;    ytlen = 0.04
;
;    plot, [0], [0], /nodata, $
;      pos=subsubplot(pos,pos3), /noerase, $
;      xrange=xrange4, xsty=1, xminor=xminor4, xtickint=xtickint4, $
;      xtickname=replicate(' ',7), xticklen=xtlen, $
;      yrange=yrange3, ysty=1, yticklen=ytlen, $
;      ytickname=replicate(' ',7), yminor=ym, ytickint=ytint
;    ;axis, yaxis=1, yrange=yrange3, ysty=1, yticklen=ytlen, $
;    ;  yminor=ym, ytickint=ytint
;    axis, xaxis=1, xrange=xrange4, xsty=1, xminor=xminor4, xticklen=xtlen, $
;      xtickint=xtint4, xtitle=xtitle5, charsize=0.8
;
;    if phot[jj].galaxy eq 'J0905+5759' then begin
;      z_bg = 0.7116d    
;      hires = rsex(path+'J0905a.dat')
;      hires_2796 = 3.d5*(hires.wave-2796.3542699d*(1.+z_bg))/$
;        (2796.3542699d*(1.+z_bg))
;      hires_2803 = 3.d5*(hires.wave-2803.5314853d*(1.+z_bg))/$
;        (2803.5314853d*(1.+z_bg))
;      oplot, hires_2796, hires.flux, psym=10, color=djs_icolor('blue'), thick=6
;      oplot, hires_2803, hires.flux, psym=10, color=djs_icolor('red'), thick=6
;
;    endif 
;
  endfor


;----------------------------------
; plot radial profile of J0905+5759
;----------------------------------

science = rsex(path+'j0905_science.txt')
star1 = rsex(path+'j0905_star1.txt')
star2 = rsex(path+'j0905_star2.txt')
star3 = rsex(path+'j0905_star3.txt')
star4 = rsex(path+'j0905_star4.txt')
star5 = rsex(path+'j0905_star5.txt')
star6 = rsex(path+'j0905_star6.txt')

avzeroprof = mean([star1[0].profile, $
                   star2[0].profile, $
                   star3[0].profile, $
                   star4[0].profile, $
                   star5[0].profile, $
                   star6[0].profile])

scale = avzeroprof / science[0].profile


fit1 = rsex(path+'j0905_sersic_neq4.txt')
fit2 = rsex(path+'j0905_sersic_neq4_reeq1.txt')
fit3 = rsex(path+'j0905_sersic_neq4_reeq2.txt')

;red, omega0=0.3, omegalambda=0.7, h100=0.7;, /verb
kpcarc = dangular(phot[jj].z, /kpc)/206265.d


  pix = 0.02
  sthick = 1
  fthick = 10

  xrange1=[0., 10.]*pix
  xrange2=xrange1*kpcarc

  plot, science.rd*pix, science.profile * scale, $
    thick=15, charsize=asize, $
    pos=subplot(32,32,[627,1023], egap=0.07), /noerase, $
    xrange=[0., 10.]*pix, xsty=5, $
    xtitle='Radius [arcsec]', $ 
    xmargin=[11., 0.], $
    yrange=[1.d-4, 5.d-2], /ylog, ysty=1, $
    ;ytitle='Surface Brightness', $
    ymargin=[4.,4.], $
    ytickname=[textoidl('10^{-4}'), textoidl('10^{-3}'), textoidl('10^{-2}')]
  axis, xaxis=0, xrange=xrange1, charsize=asize, $
    xtitle='Radius [arcsec]', xsty=1
  axis, xaxis=1, xrange=xrange2, charsize=asize, xsty=1, $
    xtitle='Radius [kpc]'
  xyouts, -0.027, 2.2d-3, 'Surface Brightness', orientation=90, charsize=asize, align=0.5

  colors = ['blue']
  colors = replicate(colors,6)
  sstyle=0

  oplot, star1.rd*pix, star1.profile, $
    color=djs_icolor(colors[0]), thick=sthick, linestyle=sstyle
  oplot, star2.rd*pix, star2.profile, $
    color=djs_icolor(colors[1]), thick=sthick, linestyle=sstyle
  oplot, star3.rd*pix, star3.profile, $
    color=djs_icolor(colors[2]), thick=sthick, linestyle=sstyle
  oplot, star4.rd*pix, star4.profile, $
    color=djs_icolor(colors[3]), thick=sthick, linestyle=sstyle
  oplot, star5.rd*pix, star5.profile, $
    color=djs_icolor(colors[4]), thick=sthick, linestyle=sstyle
  oplot, star6.rd*pix, star6.profile, $
    color=djs_icolor(colors[5]), thick=sthick, linestyle=sstyle

  oplot, fit1.rd*pix, fit1.profile * avzeroprof / fit1[0].profile, $
    color=djs_icolor('red'), thick=fthick, linestyle=2
  oplot, fit3.rd*pix, fit3.profile * avzeroprof / fit3[0].profile, $
    color=djs_icolor('dark green'), thick=fthick, linestyle=1

  csize=0.9

;legend, ['J0905+5759 profile', 'Stars in WFC3 image', 'Sersic profile (n=4, re=0.013")', 'Sersic profile (n=4, re=0.020")'], $
legend, ['J0905+5759 profile', 'Stars in WFC3 image', textoidl('Sersic profile (n=4, r_e=100 pc)'), $
  textoidl('Sersic profile (n=4, r_e=290 pc)')], $
  linestyle=[0,0,2,1], $ 
  color=[djs_icolor('black'), djs_icolor('blue'), djs_icolor('red'), djs_icolor('dark green')], $
  charsize=csize, box=0, thick=[13,5,10,10], spacing=1.5, $
  pos=[0.0, 9.d-4]  


    dfpsclose

  !p.thick=1
  !x.thick=1
  !y.thick=1
  !p.charthick=1


stop
return
end
