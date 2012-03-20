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
      symsize=1.5, color=im_color('firebrick'), $
      errcolor=im_color('firebrick'), errthick=8     


  ; put galaxy name on plot
    xyouts, 1500., 14.5, phot[jj].galaxy, size=1.5
  ; put galfit results on plot
    xyouts, 0.225, 0.82, 'data', size=1.5, align=0.5, /normal
    xyouts, 0.375, 0.82, 'model', size=1.5, align=0.5, /normal
    xyouts, 0.525, 0.82, 'residual', size=1.5, align=0.5, /normal
  ; put spectrum on plot
    if phot[jj].galaxy eq 'J0905+5759' then begin
      j0905path = j0905_path()
      lris = mrdfits(j0905path+'spectra/j0905+5759_lris_flux_v120216.fits',1)
    scale = 1D17
    yrange = [0.05,12]
    xrange1 = [4000,8700]
    specpos = [0.25, 0.15, 0.95, 0.37]

    xrange2 = xrange1/(1.0+phot[jj].z)
    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
      xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, position=specpos, $
      /noerase, color=djs_icolor('black')
    axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2
    
    djs_oplot, lris.wavelength, smooth(scale*lris.flux,1), color='grey'       
    
    plot_sigmasfr_galfit, path+'imgblock_j0905+5759_sersic_neq4.fits', $
      npix=50.
    
    endif else begin
      xyouts, 6000., 21., 'blue end', size=1.5
      xyouts, 120000., 21., 'red end', size=1.5
    endelse

    if phot[jj].galaxy eq 'J1506+5402' then begin
      plot_sigmasfr_galfit, path+'imgblock_j1506+5402_sersic_neq4.fits', $
        99., 103.
      
  
    endif

  endfor

    dfpsclose

return
end
