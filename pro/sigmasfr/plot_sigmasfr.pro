pro plot_sigmasfr
; jm12mar19ucsd - make the sigma-sfr vs mass plot
    path = getenv('HIZEA_DATA')+'/sigmasfr/'

    hst = rsex(path+'hst_sample.dat')
    kcorr = mrdfits(path+'sigmasfr_kcorrect.fits.gz',1)
    lir = mrdfits(path+'sigmasfr_lir.fits.gz',1)

    veilleux = rsex(path+'06veilleux.dat')

    xrange = [-22,-26]
    yrange = [-2,5]
    sigmasfr_edd = alog10(5000)
    
    psfile = path+'sigmasfr.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=2, $
      xmargin=[1.5,0.4], width=6.6, height=5.3
    djs_plot, [0], [0], /nodata, position=pos, xsty=5, ysty=5, $
      xrange=xrange, yrange=yrange
    polyfill, [xrange,reverse(xrange)], [xrange*0+sigmasfr_edd,xrange*0+yrange[1]], $
      /data, color=im_color('grey90'), noclip=0, /line_fill, orientation=135, spacing=0.15
    polyfill, [xrange,reverse(xrange)], [xrange*0+sigmasfr_edd,xrange*0+yrange[1]], $
      /data, color=im_color('grey90'), noclip=0, /line_fill, orientation=45, spacing=0.15
    djs_oplot, !x.crange, sigmasfr_edd*[1,1], line=0, color=im_color('grey90')
    djs_plot, [0], [0], /noerase, position=pos, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange, $
      xtickinterval=1, xtitle='Absolute H Magnitude', $
      ytitle='log (\Sigma_{SFR} / M'+sunsymbol()+' yr^{-1} kpc^{-2})'

    im_legend, 'Eddington-limited Star Formation', /left, /top, box=0, $
      charsize=1.8
    im_legend, ['HST-WISE Sample','Veilleux+06'], /left, /bottom, box=0, $
      charsize=1.5, psym=[16,15], color=['black','orange']
    
; HST sample
    good = where(lir.sfr_chary gt -900)    
    sigmasfr = alog10(lir[good].sfr_chary/(!pi*hst[good].r_e^2)/2.0)
    absmag_h = kcorr[good].k_fnuvugrizjhk_absmag_00[8]
;   niceprint, absmag_h, sigmasfr
    djs_oplot, absmag_h, sigmasfr, psym=symcat(16), symsize=1.5
    
; Veilleux+06
    sigmasfr = alog10(4.5D-44*10.0^veilleux.lir*im_lsun()/(!pi*veilleux.r50^2)/2.0)
    djs_oplot, veilleux.absmag_h, sigmasfr, psym=symcat(15), color='orange'
    
    im_plotconfig, psfile=psfile, /psclose, /pdf, /pskeep

return
end
