pro plot_seds_scuba2
;+ 
;  amd120605 -- written to make a plot showing three SEDs
;
;
;-

    path = getenv('HIZEA_DATA')+'/sigmasfr/'
    hstpath = sigmasfr_path(/multidrizzle)

    phot = mrdfits(path+'sigmasfr_photometry.fits.gz', 1)
    hst = rsex(path+'hst_sample.dat')

  readcol, path+'/swire/M82_template_norm.sed', m82wave, m82flambda
  readcol, path+'/swire/Arp220_template_norm.sed', arp220wave, arp220flambda
  readcol, path+'/swire/torus_template_norm.sed', toruswave, torusflambda

  fmt = 'f9.2,1x,f7.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.3,1x,f8.3,1x,f8.3,1x,f9.3,1x,f9.3,1x,f8.2,1x,f8.2'
  readfmt, path+'table4_avgtemplates.txt', fmt, wave, l0975, l1000, l1025, l1050, l1075, l1100, l1125, l1150, l1175, l1200, l1225, l1250, l1275, l1300, skipline=26

  restore, path+'chary_elbaz_codes/chary_elbaz.save'

; restore the iSEDfit SEDs
  isedpath = sigmasfr_path(/ised)
  isedfit_sfhgrid_dir = sigmasfr_path(/monte)
  sfhgrid_paramfile = getenv('HIZEA_DIR')+'/pro/sigmasfr/sigmasfr_sfhgrid.par'
  paramfile = isedpath+'sigmasfr_supergrid01_isedfit.par'
  model = isedfit_restore(paramfile,isedfit,iopath=isedpath,$
    isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,index=[0,1,2]) ; grab the first three

  csize=1.5

; make the plot  
    psfile = path+'seds_scuba2.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=1.5, $
      xmargin=[1.5,0.4], width=6.6, height=5.3

    filters = sigmasfr_filterlist()
    weff = k_lambda_eff(filterlist=filters) / 10000.
    isw2 = (where(strtrim(filters,2) eq 'wise_w2.par'))[0]
    isw3 = (where(strtrim(filters,2) eq 'wise_w3.par'))[0]
    isw4 = (where(strtrim(filters,2) eq 'wise_w4.par'))[0]

    xtitle = textoidl('Rest Wavelength (\mum)')
    ;ytitle = textoidl('Apparent AB Magnitude + offset')
    ytitle = textoidl('Apparent AB Magnitude')
    yrange = [23.0,8.0]
    xrange = [0.06,1000.]
    ;ticks = [0.1,0.2,0.4,0.8,1.6,3.0,6.0,12.0,24.0]
    ;xtickname=['0.1','0.2','0.4','0.8','1.6','3','6','12','24']
    ticks = [0.1,0.3,1.0,3.0,10.,30.,100.,500.]
    xtickname=['0.1','0.3','1.0','3.0','10','30','100','500']


    ;posmain = [0.12, 0.15, 0.75, 0.9]
    ;posmain = [0.12, 0.15, 0.55, 0.9]
    posmain = [0.09, 0.15, 0.51, 0.9]

    plot, [0], [0], /nodata, pos=posmain, $
      xrange=xrange, xsty=1, /xlog, $
      xtitle=xtitle, $
      xticks=n_elements(ticks)-1, xtickv=ticks, $
      ;xtickname=replicate(' ',8), $
      xtickname=xtickname, $
      yrange=yrange, ysty=1, $
      ;ytickname=replicate(' ',8), $
      ytitle=ytitle
    
  ;offset = [20.d, 0.d, 10.d]
  offset = [0., 0., 0.]

; plot seds
  for jj=0L, 0L do begin

; isedfit model
     rel = where(model[jj].wave/1D4 lt 4.3, nrel)
     djs_oplot, model[jj].wave[rel]/1D4/(1.+phot[jj].z), $
       (model[jj].flux[rel]+0.1)-offset[jj], color=im_color('grey40')
     
    ;xyouts, 50., 19.-offset[jj], phot[jj].galaxy, charsize=1.5
    ;xyouts, 0.4, 20.4-offset[jj], phot[jj].galaxy, charsize=1.8
    ;xyouts, 0.15, 22.0-offset[jj], textoidl('\Sigma_{SFR}=')+$
    ;  string(round(hst[jj].wise_sfr / (2.*!pi*hst[jj].r_e^2)/ 10.)*10., format='(i4)')+$
    ;  ' M'+sunsymbol()+textoidl(' yr^{-1} kpc^{-2}'), charsize=1.5

    if phot[jj].galaxy eq 'J0905+5759' then begin
      phot[jj].maggies[9] = phot[jj].maggies[9] * 0.8
      phot[jj].maggies[10] = phot[jj].maggies[10] * 0.8
    endif

    used = where((phot[jj].maggies gt 0.0) and $ 
      (phot[jj].ivarmaggies gt 0.0),nused)

    mab =  maggies2mag(phot[jj].maggies[used],$
      ivar=phot[jj].ivarmaggies[used],magerr=mab_err)


    oploterror, weff[used]/(1.+phot[jj].z), mab-offset[jj], mab_err, $
      psym=symcat(16), $
      symsize=1.5, color=djs_icolor('red'), $
      errcolor=djs_icolor('red'), errthick=8  



;;; plot Chary & Elbaz SED scaled to WISE ch4
;  l1mag = -2.5*alog10(nulnuinlsun[104,*] * lambda)
;  sort = sort(abs(lambda - weff[isw4]/(1.+phot[jj].z)))
;
;  loc = where(weff[used] eq weff[isw4])
;  norm = mab[loc]-offset[jj] - l1mag[sort[0]]
;  good = where(lambda gt 12.0)
;  oplot, lambda[good], l1mag[good]+norm[0], $
;   color=djs_icolor('blue'), thick=8, linestyle=2
;
;;; plot Chary & Elbaz SED scaled to WISE ch4
;  l1mag = -2.5*alog10(nulnuinlsun[50,*] * lambda)
;  sort = sort(abs(lambda - weff[isw4]/(1.+phot[jj].z)))
;
;  loc = where(weff[used] eq weff[isw4])
;  norm = mab[loc]-offset[jj] - l1mag[sort[0]]
;  good = where(lambda gt 12.0)
;  oplot, lambda[good], l1mag[good]+norm[0], $
;   color=djs_icolor('red'), thick=8, linestyle=2
;
;
;for i=0l, n_elements(lir)-1L do begin
;  l1mag = -2.5*alog10(nulnuinlsun[i,*] * lambda)
;  sort = sort(abs(lambda - weff[isw4]/(1.+phot[jj].z)))
;
;  loc = where(weff[used] eq weff[isw4])
;  norm = mab[loc]-offset[jj] - l1mag[sort[0]]
;  good = where(lambda gt 12.0)
;  oplot, lambda[good], l1mag[good]+norm[0], $
;   color=djs_icolor('red'), thick=8, linestyle=2
;
;
;endfor

;;; plot 9.75 SED scaled to WISE ch2
;  l0975mag = -2.5*alog10(l0975)
;  sort = sort(abs(wave - weff[isw3]/(1.+phot[jj].z)))
;
;  loc = where(weff[used] eq weff[isw3])
;  norm = mab[loc]-offset[jj] - l0975mag[sort[0]]
;  oplot, wave, l0975mag+norm[0], $
;   color=djs_icolor('blue'), thick=8, linestyle=2
;
;;; plot 9.75 SED scaled to WISE ch2
;  l1275mag = -2.5*alog10(l1275)
;  sort = sort(abs(wave - weff[isw3]/(1.+phot[jj].z)))
;
;  loc = where(weff[used] eq weff[isw3])
;  norm = mab[loc]-offset[jj] - l1275mag[sort[0]]
;  oplot, wave, l1275mag+norm[0], $
;   color=djs_icolor('blue'), thick=8, linestyle=2


; plot M82 SED scaled to WISE ch2
  m82mag = -2.5*alog10(m82flambda * (m82wave)^2)
  sort = sort(abs(m82wave/10000. - weff[isw2]/(1.+phot[jj].z)))

  loc = where(weff[used] eq weff[isw2])
  norm = mab[loc]-offset[jj] - m82mag[sort[0]]
  good = where(m82wave gt 43000./(1.+phot[jj].z), mgood)
  oplot, m82wave[good]/10000., m82mag[good]+norm[0], $
   color=djs_icolor('blue'), thick=8, linestyle=2
  m82norm = norm

; plot ARP220 SED scaled to WISE ch2
  arp220mag = -2.5*alog10(arp220flambda * (arp220wave)^2)
  sort = sort(abs(arp220wave/10000. - weff[isw2]/(1.+phot[jj].z)))
  loc = where(weff[used] eq weff[isw2])
  norm = mab[loc] - arp220mag[sort[0]]
  ;good = where(arp220wave gt 25000., mgood)
  good = where(arp220wave gt 43000./(1.+phot[jj].z), mgood)
  oplot, arp220wave[good]/10000., arp220mag[good]+norm[0]-offset[jj], $
    color=djs_icolor('dark green'), linestyle=3, thick=8
  arp220norm = norm

; plot QSO2 SED scaled to WISE ch2
  torusmag = -2.5*alog10(torusflambda * (toruswave)^2)
  sort = sort(abs(toruswave/10000. - weff[isw2]/(1.+phot[jj].z)))
  loc = where(weff[used] eq weff[isw2])
  norm = mab[loc] - torusmag[sort[0]]
  good = where(toruswave gt 43000./(1.+phot[jj].z), mgood)
  oplot, toruswave[good]/10000., torusmag[good]+norm[0]-offset[jj], $
    color=djs_icolor('red'), linestyle=1, thick=8

    oploterror, weff[used]/(1.+phot[jj].z), mab-offset[jj], mab_err, $
      psym=symcat(16), $
      symsize=1.5, color=djs_icolor('red'), $
      errcolor=djs_icolor('red'), errthick=8  



  endfor

  legend, ['stellar population fit', 'M82 starburst', 'Arp 220 starburst', 'obscured quasar'], $
  ;legend, ['Stellar Population Fit', 'M82 Starburst', 'Arp 220 Starburst'], $
    linestyle=[0, 2, 3, 1], $
    ;linestyle=[0, 1, 2], $
    color=[djs_icolor('black'), djs_icolor('blue'), djs_icolor('dark green'), djs_icolor('red')], $
    ;color=[djs_icolor('black'), djs_icolor('blue'), djs_icolor('dark green')], $
    charsize=1.3, box=0, pos=[xrange[0]*1.5, yrange[0]-3.2], thick=8

   lpos = 1.5

   xyouts, lpos, 11.0, 'WITHOUT SCUBA2', align=0.5
   xyouts, 70., 11.0, '?'
   
; plot second panel
   ;posher = [0.12, 0.15, 0.75, 0.9]
   ;posher = [0.55, 0.15, 0.98, 0.9]
   posher = [0.51, 0.15, 0.92, 0.9]

    plot, [0], [0], /nodata, pos=posher, /noerase, $
      xrange=xrange, xsty=1, /xlog, $
      xtitle=xtitle, $
      xticks=n_elements(ticks)-1, xtickv=ticks, $
      ;xtickname=replicate(' ',8), $
      xtickname=xtickname, $
      yrange=yrange, ysty=5, $
      ytickname=replicate(' ',8);, $
      ;ytitle=ytitle
    
  ;offset = [20.d, 0.d, 10.d]
  ;offset = [0., 0., 0.]

; plot seds
  for jj=0L, 0L do begin

; isedfit model
     rel = where(model[jj].wave/1D4 lt 4.3, nrel)
     djs_oplot, model[jj].wave[rel]/1D4/(1.+phot[jj].z), $
       model[jj].flux[rel]+0.1-offset[jj], color=im_color('grey40')
     
    ;xyouts, 50., 19.-offset[jj], phot[jj].galaxy, charsize=1.5
    ;xyouts, 0.4, 20.4-offset[jj], phot[jj].galaxy, charsize=1.8
    ;xyouts, 0.15, 22.0-offset[jj], textoidl('\Sigma_{SFR}=')+$
    ;  string(round(hst[jj].wise_sfr / (2.*!pi*hst[jj].r_e^2)/ 10.)*10., format='(i4)')+$
    ;  ' M'+sunsymbol()+textoidl(' yr^{-1} kpc^{-2}'), charsize=1.5

    if phot[jj].galaxy eq 'J0905+5759' then begin
      phot[jj].maggies[9] = phot[jj].maggies[9] * 0.8
      phot[jj].maggies[10] = phot[jj].maggies[10] * 0.8
    endif

    used = where((phot[jj].maggies gt 0.0) and $ 
      (phot[jj].ivarmaggies gt 0.0),nused)

    mab =  maggies2mag(phot[jj].maggies[used],$
      ivar=phot[jj].ivarmaggies[used],magerr=mab_err)


    oploterror, weff[used]/(1.+phot[jj].z), mab-offset[jj], mab_err, $
      psym=symcat(16), $
      symsize=1.5, color=djs_icolor('red'), $
      errcolor=djs_icolor('red'), errthick=8  


;;; plot M82 SED scaled to WISE ch2
;  m82mag = -2.5*alog10(m82flambda * (m82wave)^2)
;  sort = sort(abs(m82wave/10000. - weff[isw2]/(1.+phot[jj].z)))
;
;  loc = where(weff[used] eq weff[isw2])
;  norm = mab[loc]-offset[jj] - m82mag[sort[0]]
;  good = where(m82wave gt 43000./(1.+phot[jj].z), mgood)
;  oplot, m82wave[good]/10000., m82mag[good]+norm[0], $
;   color=djs_icolor('blue'), thick=8, linestyle=2


; plot ARP220 SED scaled to WISE ch2
  arp220mag = -2.5*alog10(arp220flambda * (arp220wave)^2)
  sort = sort(abs(arp220wave/10000. - weff[isw2]/(1.+phot[jj].z)))
  loc = where(weff[used] eq weff[isw2])
  norm = mab[loc] - arp220mag[sort[0]]
  ;good = where(arp220wave gt 25000., mgood)
  good = where(arp220wave gt 43000./(1.+phot[jj].z), mgood)
  oplot, arp220wave[good]/10000., arp220mag[good]+norm[0]-offset[jj], $
    color=im_color('grey40'), linestyle=0;, thick=6

;; plot QSO2 SED scaled to WISE ch2
;  torusmag = -2.5*alog10(torusflambda * (toruswave)^2)
;  sort = sort(abs(toruswave/10000. - weff[isw2]/(1.+phot[jj].z)))
;  loc = where(weff[used] eq weff[isw2])
;  norm = mab[loc] - torusmag[sort[0]]
;  good = where(toruswave gt 43000./(1.+phot[jj].z), mgood)
;  oplot, toruswave[good]/10000., torusmag[good]+norm[0]-offset[jj], $
;    color=djs_icolor('red'), linestyle=1, thick=8

    oploterror, weff[used]/(1.+phot[jj].z), mab-offset[jj], mab_err, $
      psym=symcat(16), $
      symsize=1.5, color=djs_icolor('red'), $
      errcolor=djs_icolor('red'), errthick=8  

    axis, yaxis=0, yrange=yrange, ysty=1, ytickname=replicate(' ',8)

    ;ytickname=[textoidl('10^{-2}'), textoidl('10^{-1}'), $
    ;  textoidl('10^{0}'), textoidl('10^{1}'), textoidl('10^{2}'), $
    ;  textoidl('10^{3}')]

    ytickname=['0.01', '0.1', '1', '10', textoidl('10^{2}'), $
      textoidl('10^{3}')]

    axis, yaxis=1, yrange=10.^(yrange/2.5*(-1.))*3631.d-23 * 1.d26, $
      /ylog, ysty=1, ytickname=ytickname;, ytitle='Flux [mJy]'

    xyouts, 0.99, 0.44, 'Flux [mJy]', orientation=90, /normal

; plot 450 micron point
    sort = sort(abs(arp220wave/10000. - 450. / (1.+phot[jj].z)))
    oplot, [0, arp220wave[sort[0]]]/10000., [0., arp220mag[sort[0]]+arp220norm[0] ], psym=6
    splog, 'Arp 220, 450: ', 10.^(-1.*(arp220mag[sort[0]]+arp220norm[0])/2.5)*3631.d3

    oploterror, [0, arp220wave[sort[0]]]/10000., [0., arp220mag[sort[0]]+arp220norm[0] ], $
      [1., 0.1], $
      psym=symcat(16), $
      symsize=1.5, color=djs_icolor('red'), $
      errcolor=djs_icolor('red'), errthick=8

; plot 850 micron point
    sort = sort(abs(arp220wave/10000. - 850. / (1.+phot[jj].z)))
    oplot, [0, arp220wave[sort[0]]]/10000., [0., arp220mag[sort[0]]+arp220norm[0] ], psym=6
    splog, 'Arp220, 850: ', 10.^(-1.*(arp220mag[sort[0]]+arp220norm[0])/2.5)*3631.d3

    oploterror, [0, arp220wave[sort[0]]]/10000., [0., arp220mag[sort[0]]+arp220norm[0] ], $
      [1., 0.1], $
      psym=symcat(16), $
      symsize=1.5, color=djs_icolor('red'), $
      errcolor=djs_icolor('red'), errthick=8

    ;loc1 = where(arp220wave 
    ;oplot, [100., 160.]/(1.+phot[jj].z), [10., 10.]

  endfor

;  legend, ['stellar population fit', 'M82 starburst', 'Arp 220 starburst', 'obscured quasar'], $
;  ;legend, ['Stellar Population Fit', 'M82 Starburst', 'Arp 220 Starburst'], $
;    linestyle=[0, 2, 3, 1], $
;    ;linestyle=[0, 1, 2], $
;    color=[djs_icolor('black'), djs_icolor('blue'), djs_icolor('dark green'), djs_icolor('red')], $
;    ;color=[djs_icolor('black'), djs_icolor('blue'), djs_icolor('dark green')], $
;    charsize=1.5, box=0, pos=[xrange[0], yrange[1]+0.5], thick=8


   xyouts, lpos, 11.0, 'WITH SCUBA2', align=0.5
   ;xyouts, 70., 11.0, '!'


;  for jj=0L, 2L do begin
;
;    if phot[jj].galaxy eq 'J0905+5759' then begin
;      hpos = [0.77, 0.15, 0.98, 0.40]
;      ;file = path+'j0905+5759_drz_sci.fits'
;      file = file_search(hstpath+'cycle1?/j0905+5759/j0905+5759_drz_sci.fits')
;      plot_sigmasfr_galfit, file, 5955., 5880., $
;        npix=120., posa=hpos, /onlydata, /ext0
;    endif
;
;    if phot[jj].galaxy eq 'J1341-0321' then begin
;      hpos = [0.77, 0.40, 0.98, 0.65]
;      ;file = path+'j1341-0321_drz_sci.fits'
;      file = file_search(hstpath+'cycle1?/j1341-0321/j1341-0321_drz_sci.fits')
;      plot_sigmasfr_galfit, file, 5258., 5457., $
;        npix=120., posa=hpos, /onlydata, /ext0
;    endif
;
;    if phot[jj].galaxy eq 'J1506+5402' then begin
;      hpos = [0.77, 0.65, 0.98, 0.90]
;      ;file = path+'j1506+5402_drz_sci.fits'
;      file = file_search(hstpath+'cycle1?/j1506+5402/j1506+5402_drz_sci.fits')
;      plot_sigmasfr_galfit, file, 5684., 5524., $
;        npix=120., posa=hpos, /onlydata, /ext0 
;    endif
;
;    xyouts, mean([hpos[0],hpos[2]]), hpos[3]-0.04, $
;      phot[jj].galaxy, charsize=1.5, /norm, align=0.5
;
;    kpcarc = dangular(phot[jj].z, /kpc)/206265.d
;
;    hst_xrange = [-120.5*kpcarc*0.02, 120.5*kpcarc*0.02]
;    hst_yrange = [-120.5*kpcarc*0.02, 120.5*kpcarc*0.02]
;    plot, [0], [0], pos=hpos, /noerase, $
;      color=djs_icolor('black'), xsty=5, ysty=5, $
;      xrange=hst_xrange, $
;      yrange=hst_yrange
;    plots, !x.crange, [hst_yrange[0], hst_yrange[0]]
;    plots, !x.crange, [hst_yrange[1], hst_yrange[1]]
;    plots, [hst_xrange[0], hst_xrange[0]], !y.crange
;    plots, [hst_xrange[1], hst_xrange[1]], !y.crange
;    plots, [5., 15.], $
;      [hst_yrange[0]*0.75, hst_yrange[0]*0.75], $
;      color=djs_icolor('black'), thick=10
;    xyouts, 10., hst_yrange[0]*0.65, '10 kpc', charsize=1.1, $
;      align=0.5, charthick=5.
;
;  endfor

    im_plotconfig, psfile=psfile, /psclose, /pdf, /pskeep

stop
return
end
