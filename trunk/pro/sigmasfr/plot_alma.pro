pro plot_alma
;+
; amd120711 -- written to make a figure for J0944+0930 for ALMA proposal
;
;-

  path = getenv('HIZEA_DATA')+'/sigmasfr/'
 
  phot = mrdfits(path+'sigmasfr_photometry.fits.gz', 1)
  hst = rsex(path+'hst_sample.dat')

  readcol, path+'/swire/M82_template_norm.sed', m82wave, m82flambda
  readcol, path+'/swire/Arp220_template_norm.sed', arp220wave, arp220flambda
  readcol, path+'/swire/torus_template_norm.sed', toruswave, torusflambda

; restore the iSEDfit SEDs
  isedpath = sigmasfr_path(/ised)
  isedfit_sfhgrid_dir = sigmasfr_path(/monte)
  sfhgrid_paramfile = getenv('HIZEA_DIR')+'/pro/sigmasfr/sigmasfr_sfhgrid.par'
  paramfile = isedpath+'sigmasfr_supergrid01_isedfit.par'
  model = isedfit_restore(paramfile,isedfit,iopath=isedpath,$
    isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,index=[0,1,2,3]) ; grab the 


; redshift for object of interest
    z = 0.608

    psfile = path+'alma.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=2, $
      xmargin=[1.5,0.4], width=6.6, height=5.3


   filters = sigmasfr_filterlist()
    weff = k_lambda_eff(filterlist=filters) / 10000.
    isw2 = (where(strtrim(filters,2) eq 'wise_w2.par'))[0]

    xtitle = textoidl('Rest Wavelength (\mum)')
    ytitle = textoidl('Apparent AB Magnitude')
    yrange = [30.0,8.0]
    xrange = [0.05,30.]
    ticks = [0.1,0.2,0.4,0.8,1.6,3.0,6.0,12.0,24.0]
    xtickname=['0.1','0.2','0.4','0.8','1.6','3','6','12','24']

    posmain = [0.12, 0.15, 0.9, 0.9]

    plot, [0], [0], /nodata, pos=posmain, /noerase, $
      xrange=xrange, xsty=1, /xlog, $
      xtitle=xtitle, $
      xticks=n_elements(ticks)-1, xtickv=ticks, $
      ;xtickname=replicate(' ',8), $
      xtickname=xtickname, $
      yrange=yrange, ysty=1, $
      ;ytickname=replicate(' ',8), $
      ytitle=ytitle

    jj = 3L

; isedfit model
     rel = where(model[jj].wave/1D4 lt 4.3, nrel)
     djs_oplot, model[jj].wave[rel]/1D4/(1.+phot[jj].z), $
       model[jj].flux[rel]+0.1;, color=im_color('grey40')


    used = where((phot[jj].maggies gt 0.0) and $ 
      (phot[jj].ivarmaggies gt 0.0),nused)

    mab =  maggies2mag(phot[jj].maggies[used],$
      ivar=phot[jj].ivarmaggies[used],magerr=mab_err)


    oploterror, weff[used]/(1.+phot[jj].z), mab, mab_err, $
      psym=symcat(16), $
      symsize=1.5, color=djs_icolor('red'), $
      errcolor=djs_icolor('red'), errthick=8  


;; plot M82 SED scaled to WISE ch2
  m82mag = -2.5*alog10(m82flambda * (m82wave)^2)
  sort = sort(abs(m82wave/10000. - weff[isw2]/(1.+phot[jj].z)))

  loc = where(weff[used] eq weff[isw2])
  norm = mab[loc] - m82mag[sort[0]]
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
  oplot, arp220wave[good]/10000., arp220mag[good]+norm[0], $
    color=djs_icolor('dark green'), linestyle=3, thick=8

; plot QSO2 SED scaled to WISE ch2
  torusmag = -2.5*alog10(torusflambda * (toruswave)^2)
  sort = sort(abs(toruswave/10000. - weff[isw2]/(1.+phot[jj].z)))
  loc = where(weff[used] eq weff[isw2])
  norm = mab[loc] - torusmag[sort[0]]
  good = where(toruswave gt 43000./(1.+phot[jj].z), mgood)
  oplot, toruswave[good]/10000., torusmag[good]+norm[0], $
    color=djs_icolor('red'), linestyle=1, thick=8

    oploterror, weff[used]/(1.+phot[jj].z), mab, mab_err, $
      psym=symcat(16), $
      symsize=1.5, color=djs_icolor('red'), $
      errcolor=djs_icolor('red'), errthick=8  

      ;hpos = [0.65, 0.15, 0.9, 0.45]
      hpos = [0.12, 0.6, 0.38, 0.9]
      file = path+'j0944+0930_drz_sci.fits';';'j1506+5402_drz_sci.fits'
      ;file = file_search(hstpath+'cycle1?/j1506+5402/j1506+5402_drz_sci.fits')
      plot_sigmasfr_galfit, file, 5943., 5995., $
        npix=120., posa=hpos, /onlydata, /ext0 

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
    plots, [0., 10.], $
      [hst_yrange[0]*0.75, hst_yrange[0]*0.75], $
      color=djs_icolor('black'), thick=10
    xyouts, 5., hst_yrange[0]*0.65, '10 kpc', charsize=1.1, $
      align=0.5, charthick=5.




     xrange = [2500., 5200.]
     ;yrange = [0.,35.]
     yrange = [0., 20.]
     xtitle = textoidl('Rest Wavelength [')+angstrom()+']'
     ytitle = textoidl('F_{\lambda} [10^{-17} erg s^{-1} cm^{-2} ')+$
       angstrom()+textoidl('^{-1}]')

    posmain = [0.12, 0.15, 0.9, 0.9]
    ;posmain = [0.12, 0.55, 0.63, 0.9]

;    plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
;      xsty=1, ysty=1, xtitle=xtitle, ytitle=ytitle, $
;      /noerase, color=djs_icolor('black'), pos=posmain, charsize=1.8
;
;; plot main spectrum
    mmt = rsex(path+'j1506+5402.txt')
;    good = where(mmt.wave lt 6860. or $
;      (mmt.wave gt 6890. and mmt.wave lt 7165.), ngood)
;    djs_oplot, mmt[good].wave/(1.+z), $
;      smooth(mmt[good].flux,1), thick=1
;    sdss_flux = mrdfits(path+'spSpec-52370-0793-244.fit', 0, hdr)
;    sdss_wave = 10.0^(findgen(sxpar(hdr, 'NAXIS1')) * sxpar(hdr, 'COEFF1') + $
;            sxpar(hdr, 'COEFF0')) 
;    good = where(sdss_wave gt 4450.*(1.+z), ngood)
;    djs_oplot, sdss_wave[good]/(1.+z), $
;      smooth(sdss_flux[good,0], 2), thick=1

;  xtitle2 = textoidl('\lambda [')+angstrom()+']'
;  yrange2 = [0., 1.49]
;
;  xtitle3 = textoidl('v [km s^{-1}]')
;  yrange3 = [0., 1.49]
;
;  xtlen = 0.06
;  csize=1.3
;
;  xrange2 = [2765., 2805.]
;  xrange3 = [-3500., 0.]
;
;;; J1506 zoom lambda
;;  ;xrange2 = [2775., 2798.]
;;  ;lpos = [0.68, 0.70, 0.84, 0.9]
;;  lpos = [0.68, 0.5, 0.84, 0.9]
;  con = 12.0
;;
;;  plot, [0], [0], /nodata, /noerase, $
;;    xrange=xrange2, xsty=1, $;xtitle=xtitle2, $
;;    ;xtickint=10., xminor=5, $
;;    xtickint=20., xminor=4, $
;;    xticklen=xtlen, $
;;    yrange=yrange2, ysty=1, $;ytitle=ytitle, $
;;    yminor=5, ytickint=0.5, yticklen=0.08, $
;;    pos=lpos, charsize=csize
;;
;  z_bg = 0.6079d    
;;  oplot, mmt.wave/(1.+z_bg), medsmooth(mmt.flux/con, 5), psym=10
;
;; J1506 zoom velocity
;  vpos = [0.68, 0.55, 0.98, 0.90]
;
;  plot, [0], [0], /nodata, /noerase, $
;    xrange=xrange3, xsty=1, $;xtitle=xtitle3, $
;    ;xtickint=1000., xminor=5, $
;    xtickint=2000., xminor=4, $
;    xticklen=xtlen, $
;    yrange=yrange3, ysty=1, yminor=5, $
;    ytickint=0.5, yticklen=0.08, ytickname=replicate(' ',3), $
;    pos=vpos, charsize=csize
;
;  mmt_2796 = 3.d5*(mmt.wave-2796.3542699d*(1.+z_bg))/$
;    (2796.3542699d*(1.+z_bg))
;  mmt_2803 = 3.d5*(mmt.wave-2803.5314853d*(1.+z_bg))/$
;   (2803.5314853d*(1.+z_bg))
;  good = where(mmt_2796 lt -1200., ngood)
;  oplot, mmt_2796[good], mmt[good].flux/con, $
;    psym=10, color=djs_icolor('blue'), thick=6
;  good = where(mmt_2803 gt -1300., ngood)
;  oplot, mmt_2803[good], mmt[good].flux/con, $
;    psym=10, color=djs_icolor('red'), thick=6
;
;  ;xyouts, mean([lpos[0], lpos[2]]), lpos[3]-0.035, 'J1506+5402', $
;  ;  align=0.5, /normal, charsize=csize
;  xyouts, mean([vpos[0], vpos[2]]), vpos[3]-0.035, textoidl('R\sim1400'), $
;    align=0.5, /normal, charsize=csize


 
; J0944
    z = 0.5138
    xrange = [-2200,0.]
    yrange=[0.,1.3]
    xpos = -1650.
    ypos = 1.0

one = rsex(path+'j0944_order85.txt')
two = rsex(path+'j0944_order84.txt')

justone = where(one.wave lt min(two.wave), njustone)
overlap = where(one.wave ge min(two.wave) and $
                one.wave le max(two.wave), noverlap)
justtwo = where(two.wave gt max(one.wave), njusttwo)

wave = [one[justone].wave, one[overlap].wave, two[justtwo].wave]
flux = [one[justone].flux, one[overlap].flux, two[justtwo].flux]
ferr = [one[justone].ferr, one[overlap].ferr, two[justtwo].ferr]

for i=0L, noverlap-1L do begin
  reltwo = where(two.wave eq one[overlap[i]].wave)
  values = [one[overlap[i]].flux, two[reltwo].flux]
  sigmas = [one[overlap[i]].ferr, two[reltwo].ferr]
  result = amd_wav(values,sigmas, /silent)
  flux[overlap[i]] = result.value
  ferr[overlap[i]] = result.sigma
endfor

    ;spec = rsex('j0944_order85.txt')
    vel_2796 = 3.d5*(wave-2796.351d*(1.+z))/(2796.351d*(1.+z))
    vel_2803 = 3.d5*(wave-2803.531d*(1.+z))/(2796.351d*(1.+z))

    xtitle = textoidl('Outflow Velocity [km s^{-1}]')
    ytitle = 'Flux'
    ;yticklen = 0.05
    xticklen = 0.05

    plot, [0], [0], /nodata, /noerase, $
      xrange=xrange, xsty=1, xtitle=xtitle, xticklen=xticklen, $
      yticklen=yticklen, ytickinterval=0.5, yminor=5., $
      yrange=yrange, ysty=1, ytitle=ytitle, $
      ;pos=[locx[2],locy[2],locx[3],locy[3]], $
      ;pos=[0.2, 0.6, 0.68, 0.88], $
      pos = [0.35, 0.27, 0.87, 0.47], $
      charsize=1.5, thick=2.


   ; plot, vel_2796, flux, /noerase, $
   ;   xrange=xrange, xsty=5, yticklen=yticklen, $
   ;   yrange=yrange, ysty=1, ytickinterval=0.5, yminor=5., $
   ;   ;pos=[locx[2],locy[2],locx[3],locy[3]], $
   ;   pos=[0.3, 0.2, 0.8, 0.4], $
   ;   charsize=1.5, thick=2    
    good = where(vel_2796 lt -1050)
    oplot, vel_2796[good], flux[good], thick=2;, color=djs_icolor('blue')
    good = where(vel_2803 gt -1150)
    oplot, vel_2803[good], flux[good], thick=2;, color=djs_icolor('red')
    ;axis, xaxis=1, xrange=xrange, $
    ;  xsty=1, xminor=10, xtickinterval=1000, xticklen=xticklen, $
    ;  xcharsize=1.0, $
    ;  xtickname=['-2000', '-1000', ' ']

    ;axis, xaxis=0, xrange=xrange-768., $
    ;  xsty=1, xminor=10, xtickinterval=1000, xticklen=xticklen, $
    ;  xcharsize=1.5, $
    ;  xtickname=['-2000','-1000']

    ;spec = rsex('j0944_order84.txt')
    ;vel_2796 = 3.d5*(spec.wave-2796.351d*(1.+z))/(2796.351d*(1.+z))
    ;oplot, vel_2796, spec.flux, color=djs_icolor('red'), thick=1


legend, ['stellar population fit', 'M82 starburst', 'Arp 220 starburst', 'obscured quasar'], $
  ;legend, ['Stellar Population Fit', 'M82 Starburst', 'Arp 220 Starburst'], $
    linestyle=[0, 2, 3, 1], $
    ;linestyle=[0, 1, 2], $
    color=[djs_icolor('black'), djs_icolor('blue'), djs_icolor('dark green'), djs_icolor('red')], $
    ;color=[djs_icolor('black'), djs_icolor('blue'), djs_icolor('dark green')], $
    charsize=1.5, box=0, pos=[xrange[0]+150., yrange[1]+2.6], thick=8


    im_plotconfig, psfile=psfile, /psclose, /pdf, /pskeep


stop

end
