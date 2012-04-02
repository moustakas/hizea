pro plotseds_sigmasfr
; ad12mar19ucsd - plot up the SEDs

  !p.thick=4
  !x.thick=4
  !y.thick=4
  !p.charthick=4

  path = getenv('HIZEA_DATA')+'/sigmasfr/'
  phot = mrdfits(path+'sigmasfr_photometry.fits.gz', 1)
  filters = sigmasfr_filterlist()
  weff = k_lambda_eff(filterlist=filters)

  dfpsplot, path+'hst_sample_seds.ps', /landscape, /color

  for jj=0L, n_elements(phot)-1L do begin

    used = where((phot[jj].maggies gt 0.0) and $ 
      (phot[jj].ivarmaggies gt 0.0),nused)

    mab =  maggies2mag(phot[jj].maggies[used],$
      ivar=phot[jj].ivarmaggies[used],magerr=mab_err)

    xtitle1 = textoidl('Observed-Frame Wavelength (\AA)')
    xtitle2 = textoidl('Rest-Frame Wavelength (\AA)')
    ytitle1 = textoidl('Apparent AB Magnitude')
    ytitle2 = textoidl('F_{\lambda}')

    yrange = [24.0,13.5]
    fac = 1.
    ;fac = 10000.
    xrange1 = [1200,300000]/fac

    ticks1 = [2000,4000,8000,16000,30000,60000,120000,240000]/fac
    ticks2 = [1000,2000,4000,8000,16000,30000,60000,120000]/fac

    xrange2 = xrange1/(1.0+phot[jj].z)
    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
      xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, xlog=1, $ 
      position=pos, xtickformat='(I0)', xticks=n_elements(ticks1)-1, $
      xtickv=ticks1, color=djs_icolor('black'), charsize=1.5, $
      xmargin=[7.,2.], ymargin=[4.,3.]
    axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2, xlog=1, $
      xticks=n_elements(ticks2)-1, xtickv=ticks2, xtickformat='(I0)', $
      charsize=1.5

    oploterror, weff[used]/fac, mab, mab_err, psym=symcat(16), $
      symsize=1.5, color=djs_icolor('red'), $
      errcolor=djs_icolor('red'), errthick=8     

;if max(weff[used]) gt 20000. then begin

  readcol, path+'/swire/M82_template_norm.sed', m82wave, m82flambda
  m82mag = -2.5*alog10(m82flambda * (m82wave)^2)
  sort = sort(abs(m82wave - weff[8])/(1.+phot[jj].z))
  loc = where(weff[used] eq weff[8])
  offset = mab[loc] - m82mag[sort[0]]
  good = where(m82wave gt 25000., mgood)
  oplot, m82wave[good]*(1.+phot[jj].z), m82mag[good]+offset[0], color=djs_icolor('blue')
;stop
;endif

; put galaxy name on plot
    xyouts, 0.22, 0.55, phot[jj].galaxy, size=2, align=0.5, /normal, charthick=4
; put galfit results on plot
  ; horizonal spacing on the top left
    posa = [0.12, 0.65, 0.32, 0.90]
    posb = [0.32, 0.65, 0.52, 0.90] 
    posc = [0.52, 0.65, 0.72, 0.90]

    if phot[jj].galaxy eq 'J0905+5759' then $
      plot_sigmasfr_galfit, path+'imgblock_j0905+5759_sersic_neq4.fits', $
        npix=60., posa=posa, posb=posb, posc=posc
    if phot[jj].galaxy eq 'J1341-0321' then $
      plot_sigmasfr_galfit, path+'imgblock_j1341-0321_sersic_neq4.fits', $
        105., 106., posa=posa, posb=posb, posc=posc  
    if phot[jj].galaxy eq 'J1506+5402' then $
      plot_sigmasfr_galfit, path+'imgblock_j1506+5402_sersic_neq4.fits', $
        99., 103., posa=posa, posb=posb, posc=posc

    plot, [0], [0], pos=posa, /noerase, $
      color=djs_icolor('black'), $
      xtickname=replicate(' ',6), $
      ytickname=replicate(' ',6)

    plot, [0], [0], pos=posb, /noerase, $
      color=djs_icolor('black'), $
      xtickname=replicate(' ',6), $
      ytickname=replicate(' ',6)

    plot, [0], [0], pos=posc, /noerase, $
      color=djs_icolor('black'), $
      xtickname=replicate(' ',6), $
      ytickname=replicate(' ',6)

    dy = 0.03
    xyouts, mean([posa[0],posa[2]]), posa[1]-dy, 'HST WFC3/F814W', $
      size=1.5, align=0.5, /normal
    xyouts, mean([posb[0],posb[2]]), posb[1]-dy, 'GALFIT model', $
      size=1.5, align=0.5, /normal
    xyouts, mean([posc[0],posc[2]]), posc[1]-dy, 'residual', $
      size=1.5, align=0.5, /normal

; put lores spectrum on plot
  ; spectra across the bottom
    pos1 = [0.25, 0.13, 0.65, 0.35]
    pos2 = [0.65, 0.13, 0.80, 0.35] 
    pos3 = [0.80, 0.13, 0.95, 0.35]

    xrange1 = [4000,8700]
    xrange2 = xrange1/(1.0+phot[jj].z)
    yrange = [0.05,12]
    xtitle1 = textoidl('Observed-Frame Wavelength (\AA)')
    xtitle2 = textoidl('Rest-Frame Wavelength (\AA)')
    ytitle1 = textoidl('F_{\lambda}')
    xtlen=0.04

    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
      xsty=5, ysty=1, xtitle=xtitle1, ytitle=ytitle1, position=pos1, $
      /noerase, color=djs_icolor('black'), xticklen=xtlen
    axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2, xticklen=xtlen
    axis, xaxis=0, xsty=1, xrange=xrange2, xtickname=replicate(' ',6), $
        xticklen=xtlen    

    if phot[jj].galaxy eq 'J0905+5759' then begin
      j0905path = j0905_path()
      lris = mrdfits(j0905path+'spectra/j0905+5759_lris_flux_v120216.fits',1)
      scale = 1D17
      djs_oplot, lris.wavelength, smooth(scale*lris.flux,1), color='grey'       
    endif

  ; plot hires spectrum    
    xrange3 = [2765., 2785.]
    yrange3 = [0., 1.2]
    xminor3 = 4.
    xtint3 = 10.
    xtitle3 = textoidl('Rest \lambda [')+angstrom()+']'
    xtlen=0.04
    yminor3 = 5

    plot, [0], [0], /nodata, pos=pos2, /noerase, $
        xrange=xrange3, xsty=1, xminor=xminor3, xtickint=xtint3, $
        xtickname=replicate(' ',7), xticklen=xtlen, $
        ytickname=replicate(' ',7), $
        yrange=yrange3, ysty=5, yminor=yminor3, ytickint=ytint
      axis, yaxis=1, yrange=yrange3, ysty=1, yminor=yminor3, ytickint=ytint, $
        ytickname=replicate(' ',7), yticklen=ytlen
      axis, xaxis=1, xrange=xrange3, xsty=1, xminor=xminor3, $
        xtickint=xtint3, xticklen=xtlen, xtitle=xtitle3
      
    if phot[jj].galaxy eq 'J0905+5759' then begin
      z_bg = 0.7116d    
      hires = rsex(path+'J0905a.dat')
      oplot, hires.wave/(1.+z_bg), medsmooth(hires.flux, 5), psym=10
    endif 

  ; plot hires spectrum in velocity space
    xrange4 = [-3000., -2200.]
    xminor4 = 4.
    xtint4 = 400.
    xtitle4 = textoidl('Mg II Velocity [km s^{-1}]')
    ytint = 0.5
    ytm = 5
    xtlen = 0.04
    ytlen = 0.04

    plot, [0], [0], /nodata, $
      pos=pos3, /noerase, $
      xrange=xrange4, xsty=1, xminor=xminor4, xtickint=xtickint4, $
      xtickname=replicate(' ',7), xticklen=xtlen, $
      yrange=yrange3, ysty=1, yticklen=ytlen, $
      ytickname=replicate(' ',7), yminor=ym, ytickint=ytint
    ;axis, yaxis=1, yrange=yrange3, ysty=1, yticklen=ytlen, $
    ;  yminor=ym, ytickint=ytint
    axis, xaxis=1, xrange=xrange4, xsty=1, xminor=xminor4, xticklen=xtlen, $
      xtickint=xtint4, xtitle=xtitle4

    if phot[jj].galaxy eq 'J0905+5759' then begin
      z_bg = 0.7116d    
      hires = rsex(path+'J0905a.dat')
      hires_2796 = 3.d5*(hires.wave-2796.3542699d*(1.+z_bg))/$
        (2796.3542699d*(1.+z_bg))
      hires_2803 = 3.d5*(hires.wave-2803.5314853d*(1.+z_bg))/$
        (2803.5314853d*(1.+z_bg))
      oplot, hires_2796, hires.flux, psym=10, color=djs_icolor('blue'), thick=6
      oplot, hires_2803, hires.flux, psym=10, color=djs_icolor('red'), thick=6

    endif 

  endfor

    dfpsclose

  !p.thick=1
  !x.thick=1
  !y.thick=1
  !p.charthick=1

return
end
