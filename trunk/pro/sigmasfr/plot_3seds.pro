pro plot_3seds
;+ 
;  amd120406 -- written to make a plot showing three SEDs
;
;
;-

    path = getenv('HIZEA_DATA')+'/sigmasfr/'
    hstpath = sigmasfr_path(/multidrizzle)

    phot = mrdfits(path+'sigmasfr_photometry.fits.gz', 1)
    hst = rsex(path+'hst_sample.dat')

  readcol, path+'/swire/M82_template_norm.sed', m82wave, m82flambda
  readcol, path+'/swire/Arp220_template_norm.sed', arp220wave, arp220flambda
  readcol, path+'/swire/Torus_template_norm.sed', toruswave, torusflambda

; restore the iSEDfit SEDs
  isedpath = sigmasfr_path(/ised)
  isedfit_sfhgrid_dir = sigmasfr_path(/monte)
  sfhgrid_paramfile = getenv('HIZEA_DIR')+'/pro/sigmasfr/sigmasfr_sfhgrid.par'
  paramfile = isedpath+'sigmasfr_supergrid01_isedfit.par'
  model = isedfit_restore(paramfile,isedfit,iopath=isedpath,$
    isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,index=[0,1,2]) ; grab the first three

; make the plot  
    psfile = path+'3seds.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=2, $
      xmargin=[1.5,0.4], width=6.6, height=5.3

    filters = sigmasfr_filterlist()
    weff = k_lambda_eff(filterlist=filters) / 10000.
    isw2 = (where(strtrim(filters,2) eq 'wise_w2.par'))[0]

    xtitle = textoidl('Rest Wavelength (\mum)')
    ytitle = textoidl('Apparent AB Magnitude + offset')
    yrange = [25.0,1.0]
    xrange = [0.05,30.]
    ticks = [0.1,0.2,0.4,0.8,1.6,3.0,6.0,12.0,24.0]
    xtickname=['0.1','0.2','0.4','0.8','1.6','3','6','12','24']

    posmain = [0.12, 0.15, 0.75, 0.9]

    plot, [0], [0], /nodata, pos=posmain, $
      xrange=xrange, xsty=1, /xlog, $
      xtitle=xtitle, $
      xticks=n_elements(ticks)-1, xtickv=ticks, $
      ;xtickname=replicate(' ',8), $
      xtickname=xtickname, $
      yrange=yrange, ysty=1, $
      ;ytickname=replicate(' ',8), $
      ytitle=ytitle
    
  offset = [10.d, 0.d, 5.d]
  ;offset = [0., 0., 0.]

; plot seds
  for jj=0L, 2L do begin

; isedfit model
     rel = where(model[jj].wave/1D4 lt 4.3, nrel)
     djs_oplot, model[jj].wave[rel]/1D4/(1.+phot[jj].z), $
       model[jj].flux[rel]-offset[jj], color=im_color('grey40')
     
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


;; plot M82 SED scaled to WISE ch2
  m82mag = -2.5*alog10(m82flambda * (m82wave)^2)
  sort = sort(abs(m82wave/10000. - weff[isw2]/(1.+phot[jj].z)))

  loc = where(weff[used] eq weff[isw2])
  norm = mab[loc]-offset[jj] - m82mag[sort[0]]
  good = where(m82wave gt 43000./(1.+phot[jj].z), mgood)
  oplot, m82wave[good]/10000., m82mag[good]+norm[0], $
   color=djs_icolor('blue'), thick=8, linestyle=2


; plot ARP220 SED scaled to WISE ch2
  arp220mag = -2.5*alog10(arp220flambda * (arp220wave)^2)
  sort = sort(abs(arp220wave/10000. - weff[isw2]/(1.+phot[jj].z)))
  loc = where(weff[used] eq weff[isw2])
  norm = mab[loc] - arp220mag[sort[0]]
  ;good = where(arp220wave gt 25000., mgood)
  good = where(arp220wave gt 43000./(1.+phot[jj].z), mgood)
  oplot, arp220wave[good]/10000., arp220mag[good]+norm[0]-offset[jj], $
    color=djs_icolor('dark green'), linestyle=3, thick=8

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
    charsize=1.5, box=0, pos=[xrange[0], yrange[1]+0.5], thick=8, pspacing=2.5


   

  for jj=0L, 2L do begin

    if phot[jj].galaxy eq 'J0905+5759' then begin
      hpos = [0.77, 0.15, 0.98, 0.40]
      ;file = path+'j0905+5759_drz_sci.fits'
      file = file_search(hstpath+'cycle1?/j0905+5759/j0905+5759_drz_sci.fits')
      plot_sigmasfr_galfit, file, 5955., 5880., $
        npix=120., posa=hpos, /onlydata, /ext0
    endif

    if phot[jj].galaxy eq 'J1341-0321' then begin
      hpos = [0.77, 0.40, 0.98, 0.65]
      ;file = path+'j1341-0321_drz_sci.fits'
      file = file_search(hstpath+'cycle1?/j1341-0321/j1341-0321_drz_sci.fits')
      plot_sigmasfr_galfit, file, 5258., 5457., $
        npix=120., posa=hpos, /onlydata, /ext0
    endif

    if phot[jj].galaxy eq 'J1506+5402' then begin
      hpos = [0.77, 0.65, 0.98, 0.90]
      ;file = path+'j1506+5402_drz_sci.fits'
      file = file_search(hstpath+'cycle1?/j1506+5402/j1506+5402_drz_sci.fits')
      plot_sigmasfr_galfit, file, 5684., 5524., $
        npix=120., posa=hpos, /onlydata, /ext0 
    endif

    xyouts, mean([hpos[0],hpos[2]]), hpos[3]-0.04, $
      phot[jj].galaxy, charsize=1.5, /norm, align=0.5

    kpcarc = dangular(phot[jj].z, /kpc)/206265.d

    hst_xrange = [-120.5*kpcarc*0.02, 120.5*kpcarc*0.02]
    hst_yrange = [-120.5*kpcarc*0.02, 120.5*kpcarc*0.02]
    plot, [0], [0], pos=hpos, /noerase, $
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

  endfor

    im_plotconfig, psfile=psfile, /psclose, /pdf, /pskeep

; plot escape velocity v. outflow velocity

plot, 

stop
return
end
