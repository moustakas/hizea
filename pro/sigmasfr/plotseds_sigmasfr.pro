pro plotseds_sigmasfr
; ad12mar19ucsd - plot up the SEDs

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

    yrange = [24.0,13.5]
    xrange1 = [1200,300000]

    ticks1 = [2000,4000,8000,16000,30000,60000,120000,240000]
    ticks2 = [1000,2000,4000,8000,16000,30000,60000,120000]

    xrange2 = xrange1/(1.0+phot[jj].z)
    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
      xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, xlog=1, $ 
      position=pos, xtickformat='(I0)', xticks=n_elements(ticks1)-1, $
      xtickv=ticks1, color=djs_icolor('black')
    axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2, xlog=1, $
      xticks=n_elements(ticks2)-1, xtickv=ticks2, xtickformat='(I0)'

    oploterror, weff[used], mab, mab_err, psym=symcat(16), $
      symsize=1.5, color=djs_icolor('red'), $
      errcolor=djs_icolor('red'), errthick=8     


; put galaxy name on plot
    xyouts, 0.22, 0.55, phot[jj].galaxy, size=2, align=0.5, /normal, charthick=4
; put galfit results on plot
  ; horizonal spacing on the top left
    posa = [0.12, 0.67, 0.32, 0.92]
    posb = [0.32, 0.67, 0.52, 0.92] 
    posc = [0.52, 0.67, 0.72, 0.92]

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

; put spectrum on plot
    if phot[jj].galaxy eq 'J0905+5759' then begin
      j0905path = j0905_path()
      lris = mrdfits(j0905path+'spectra/j0905+5759_lris_flux_v120216.fits',1)
      scale = 1D17
      yrange = [0.05,12]
      xrange1 = [4000,8700]
    ; across the bottom 
      specpos = [0.25, 0.13, 0.65, 0.35]
      hirespos1 = [0.65, 0.13, 0.8, 0.35] 
      hirespos2 = [0.8, 0.13, 0.95, 0.35]      


      xrange2 = xrange1/(1.0+phot[jj].z)
      plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
        xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, position=specpos, $
        /noerase, color=djs_icolor('black')
      axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2
    
      djs_oplot, lris.wavelength, smooth(scale*lris.flux,1), color='grey'       
          
      z_bg = 0.7116d    
      hires = rsex(path+'J0905a.dat')
      hires_2796 = 3.d5*(hires.wave-2796.3542699d*(1.+z_bg))/$
        (2796.3542699d*(1.+z_bg))
      hires_2803 = 3.d5*(hires.wave-2803.5314853d*(1.+z_bg))/$
        (2803.5314853d*(1.+z_bg))
      
      xrange3 = [4730., 4770.]
      yrange3 = [0., 1.2]
      xminor3 = 4.
      xtickint3 = 20.
      xrange4 = [-3000., -2200.]
      xminor4 = 4.
      xtickint4 = 400.

      plot, hires.wave, medsmooth(hires.flux, 5), psym=10, thick=6, $
        pos=hirespos1, /noerase, $
        xrange=xrange3, xsty=1, xminor=xminor3, xtickint=xtickint3, $
        ytickname=replicate(' ',7), $
        yrange=yrange3, ysty=1

      plot, hires_2796, hires.flux, psym=10, thick=6, $
        pos=hirespos2, /noerase, $
        xrange=xrange4, xsty=1, xminor=xminor4, xtickint=xtickint4, $
        yrange=yrange3, ysty=1, $
        ytickname=replicate(' ',7)

      oplot, hires_2796, hires.flux, psym=10, color=djs_icolor('blue'), thick=6
      oplot, hires_2803, hires.flux, psym=10, color=djs_icolor('red'), thick=6

    endif else begin
      xyouts, 6000., 21., 'blue end', size=1.5
      xyouts, 120000., 21., 'red end', size=1.5
    endelse


      



  endfor

    dfpsclose

return
end
