pro plot_spectra
;+ 
;  amd120405 -- written to make a plot of the optical spectra
;
;
;-

    path = getenv('HIZEA_DATA')+'/sigmasfr/'
    phot = mrdfits(path+'sigmasfr_photometry.fits.gz', 1)


    psfile = path+'spectra.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=2, $
      xmargin=[1.5,0.4], width=6.6, height=5.3

     xrange = [2500., 5200.]
     yrange = [0.,35.]
     xtitle = textoidl('Rest Wavelength [')+angstrom()+']'
     ytitle = textoidl('F_{\lambda} [10^{-17} erg s^{-1} cm^{-2} ')+$
       angstrom()+textoidl('^{-1}] + offset')

    posmain = [0.12, 0.15, 0.63, 0.9]

    plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
      xsty=1, ysty=1, xtitle=xtitle, ytitle=ytitle, $
      /noerase, color=djs_icolor('black'), pos=posmain, charsize=1.8


  offset = [20., 0., 10.]
  corr = [11., 7., 10.]

; plot lores spectra
  for jj=0L, 2L do begin

    xyouts, 3800., offset[jj]+corr[jj], phot[jj].galaxy, charsize=1.8

  if phot[jj].galaxy eq 'J0905+5759' then begin
    ;offset = 0.
    j0905path = j0905_path()
    lris = mrdfits(j0905path+'spectra/j0905+5759_lris_flux_v120216.fits',1)
    scale = 1D17
    good = where(lris.wavelength lt 6860. or $
      (lris.wavelength gt 6890. and lris.wavelength lt 7590.) or $
      lris.wavelength gt 7675., ngood)
    djs_oplot, lris[good].wavelength/(1.+phot[jj].z), $
      smooth(scale*lris[good].flux,1)+offset[jj], thick=1
    sdss_flux = mrdfits(path+'spSpec-51924-0483-628.fit', 0, hdr)
    sdss_wave = 10.0^(findgen(sxpar(hdr, 'NAXIS1')) * sxpar(hdr, 'COEFF1') + $
            sxpar(hdr, 'COEFF0')) 
    good = where(sdss_wave gt 5200.*(1.+phot[jj].z), ngood)
    djs_oplot, sdss_wave[good]/(1.+phot[jj].z), $
      smooth(sdss_flux[good,0], 2)+offset[jj], thick=1
  endif

  if phot[jj].galaxy eq 'J1506+5402' then begin
    ;offset = 20.
    mmt = rsex(path+'j1506+5402.txt')
    good = where(mmt.wave lt 6860. or $
      (mmt.wave gt 6890. and mmt.wave lt 7165.), ngood)
    ;good = where(mmt.wave gt 0., ngood)
    djs_oplot, mmt[good].wave/(1.+phot[jj].z), $
      smooth(mmt[good].flux,1)+offset[jj], thick=1
    sdss_flux = mrdfits(path+'spSpec-52370-0793-244.fit', 0, hdr)
    sdss_wave = 10.0^(findgen(sxpar(hdr, 'NAXIS1')) * sxpar(hdr, 'COEFF1') + $
            sxpar(hdr, 'COEFF0')) 
    good = where(sdss_wave gt 4450.*(1.+phot[jj].z), ngood)
    djs_oplot, sdss_wave[good]/(1.+phot[jj].z), $
      smooth(sdss_flux[good,0], 2)+offset[jj], thick=1
  endif
 
  if phot[jj].galaxy eq 'J1341-0321' then begin
    ;offset=10.
    mage = rsex(path+'j1341-0321.txt')
    good = where(mage.wave lt 6860. or $
      (mage.wave gt 6890. and mage.wave lt 7595.) or $
      mage.wave gt 7670., ngood)
    djs_oplot, mage[good].wave/(1.+phot[jj].z), $
      smooth(mage[good].flux,5)+offset[jj], thick=1
    sdss_flux = mrdfits(path+'spSpec-52433-0913-300.fit', 0, hdr)
    sdss_wave = 10.0^(findgen(sxpar(hdr, 'NAXIS1')) * sxpar(hdr, 'COEFF1') + $
            sxpar(hdr, 'COEFF0')) 
    good = where(sdss_wave gt 5200.*(1.+phot[jj].z), ngood)
    djs_oplot, sdss_wave[good]/(1.+phot[jj].z), $
      smooth(sdss_flux[good,0], 2)+offset[jj], thick=1
  endif

  endfor

  xtitle2 = textoidl('\lambda [')+angstrom()+']'
  yrange2 = [0., 1.49]

  xtitle3 = textoidl('v [km s^{-1}]')
  yrange3 = [0., 1.49]

  xtlen = 0.06
  csize=1.3

  xrange2 = [2765., 2805.]
  xrange3 = [-3500., 0.]

; J0905 zoom lambda
  ;xrange2 = [2765., 2785.]
  lpos = [0.68, 0.15, 0.84, 0.35]

  plot, [0], [0], /nodata, /noerase, $
    xrange=xrange2, xsty=1, xtitle=xtitle2, $
    ;xtickint=10., xminor=5, $
    xtickint=20., xminor=4, $
    xticklen=xtlen, $
    yrange=yrange2, ysty=1, $;ytitle=ytitle, $
    yminor=5, ytickint=0.5, yticklen=0.08, $
    pos=lpos, charsize=csize

  z_bg = 0.7116d    
  hires = rsex(path+'J0905a.dat')
  oplot, hires.wave/(1.+z_bg), medsmooth(hires.flux, 7), psym=10

; J0905 zoom velocity
  ;xrange3 = [-2999., -2250.]
  vpos = [0.84, 0.15, 0.98, 0.35]
  ;ticks = [-2400., -2900.]

  plot, [0], [0], /nodata, /noerase, $
    xrange=xrange3, xsty=1, xtitle=xtitle3, $
    xticks=n_elements(ticks)-1, xtickv=ticks, $
    ;xtickint=400., xminor=4, $
    xtickint=2000., xminor=4, $
    xticklen=xtlen, $
    yrange=yrange3, ysty=1, yminor=5, $
    ytickint=0.5, yticklen=0.08, ytickname=replicate(' ',3), $
    pos=vpos, charsize=csize

  hires_2796 = 3.d5*(hires.wave-2796.3542699d*(1.+z_bg))/$
    (2796.3542699d*(1.+z_bg))
  hires_2803 = 3.d5*(hires.wave-2803.5314853d*(1.+z_bg))/$
   (2803.5314853d*(1.+z_bg))
  good = where(hires_2796 lt -2200., ngood)
  oplot, hires_2796[good], smooth(hires[good].flux, 4), psym=10, color=djs_icolor('blue'), thick=6
  good = where(hires_2803 gt -3000., ngood)
  oplot, hires_2803[good], smooth(hires[good].flux, 4), psym=10, color=djs_icolor('red'), thick=6

  ;xyouts, mean([lpos[0], lpos[2]]), lpos[3]-0.035, 'J0905+5759', $
  ;  align=0.5, /normal, charsize=csize
  xyouts, mean([vpos[0], vpos[2]]), vpos[3]-0.035, textoidl('R\sim37,000'), $
    align=0.5, /normal, charsize=csize
  ;xyouts, mean([vpos[0], vpos[2]]), vpos[3]-0.035, textoidl('FWHM\sim8 km/s'), $
  ;  align=0.5, /normal, charsize=csize
  ;xyouts, mean([vpos[0], vpos[2]]), vpos[3]-0.08, textoidl('\sim8 km/s'), $
  ;  align=0.5, /normal, charsize=csize

; J1341 zoom lambda
  ;xrange2 = [2782., 2805.]
  lpos = [0.68, 0.425, 0.84, 0.625]
  con = 13.7

  plot, [0], [0], /nodata, /noerase, $
    xrange=xrange2, xsty=1, $;xtitle=xtitle2, $
    ;xtickint=10., xminor=5, $
    xtickint=20., xminor=4, $
    xticklen=xtlen, $
    yrange=yrange2, ysty=1, $;ytitle=ytitle, $
    yminor=5, ytickint=0.5, yticklen=0.08, $
    pos=lpos, charsize=csize

  z_bg = 0.6612d    
  oplot, mage.wave/(1.+z_bg), medsmooth(mage.flux/con, 2), psym=10

; J1341 zoom velocity
  ;xrange3 = [-1400., 0.]
  vpos = [0.84, 0.425, 0.98, 0.625]

  plot, [0], [0], /nodata, /noerase, $
    xrange=xrange3, xsty=1, $;xtitle=xtitle3, $
    ;xtickint=1000., xminor=5, $
    xtickint=2000., xminor=4, $
    xticklen=xtlen, $
    yrange=yrange3, ysty=1, yminor=5, $
    ytickint=0.5, yticklen=0.08, ytickname=replicate(' ',3), $
    pos=vpos, charsize=csize

  mage_2796 = 3.d5*(mage.wave-2796.3542699d*(1.+z_bg))/$
    (2796.3542699d*(1.+z_bg))
  mage_2803 = 3.d5*(mage.wave-2803.5314853d*(1.+z_bg))/$
   (2803.5314853d*(1.+z_bg))
  good = where(mage_2796 lt -500., ngood)
  oplot, mage_2796[good], mage[good].flux/con, $
    psym=10, color=djs_icolor('blue'), thick=6
  good = where(mage_2803 gt -1000., ngood)
  oplot, mage_2803[good], mage[good].flux/con, $
    psym=10, color=djs_icolor('red'), thick=6

  ;xyouts, mean([lpos[0], lpos[2]]), lpos[3]-0.035, 'J1341-0321', $
  ;  align=0.5, /normal, charsize=csize
  xyouts, mean([vpos[0], vpos[2]]), vpos[3]-0.035, textoidl('R\sim4100'), $
    align=0.5, /normal, charsize=csize

  mage_4861 = 3.d5*(mage.wave-4861.d*(1.+z_bg))/$
    (4861.d*(1.+z_bg))
  mage_3727 = 3.d5*(mage.wave-3727.d*(1.+z_bg))/$
    (3727.d*(1.+z_bg))

; J1506 zoom lambda
  ;xrange2 = [2775., 2798.]
  lpos = [0.68, 0.70, 0.84, 0.9]
  con = 12.0

  plot, [0], [0], /nodata, /noerase, $
    xrange=xrange2, xsty=1, $;xtitle=xtitle2, $
    ;xtickint=10., xminor=5, $
    xtickint=20., xminor=4, $
    xticklen=xtlen, $
    yrange=yrange2, ysty=1, $;ytitle=ytitle, $
    yminor=5, ytickint=0.5, yticklen=0.08, $
    pos=lpos, charsize=csize

  z_bg = 0.6079d    
  oplot, mmt.wave/(1.+z_bg), medsmooth(mmt.flux/con, 5), psym=10

; J1506 zoom velocity
  ;xrange3 = [-2200., -800.]
  vpos = [0.84, 0.70, 0.98, 0.90]

  plot, [0], [0], /nodata, /noerase, $
    xrange=xrange3, xsty=1, $;xtitle=xtitle3, $
    ;xtickint=1000., xminor=5, $
    xtickint=2000., xminor=4, $
    xticklen=xtlen, $
    yrange=yrange3, ysty=1, yminor=5, $
    ytickint=0.5, yticklen=0.08, ytickname=replicate(' ',3), $
    pos=vpos, charsize=csize

  mmt_2796 = 3.d5*(mmt.wave-2796.3542699d*(1.+z_bg))/$
    (2796.3542699d*(1.+z_bg))
  mmt_2803 = 3.d5*(mmt.wave-2803.5314853d*(1.+z_bg))/$
   (2803.5314853d*(1.+z_bg))
  good = where(mmt_2796 lt -1200., ngood)
  oplot, mmt_2796[good], mmt[good].flux/con, $
    psym=10, color=djs_icolor('blue'), thick=6
  good = where(mmt_2803 gt -1300., ngood)
  oplot, mmt_2803[good], mmt[good].flux/con, $
    psym=10, color=djs_icolor('red'), thick=6

  ;xyouts, mean([lpos[0], lpos[2]]), lpos[3]-0.035, 'J1506+5402', $
  ;  align=0.5, /normal, charsize=csize
  xyouts, mean([vpos[0], vpos[2]]), vpos[3]-0.035, textoidl('R\sim1400'), $
    align=0.5, /normal, charsize=csize

    im_plotconfig, psfile=psfile, /psclose, /pdf, /pskeep




stop
return
end
