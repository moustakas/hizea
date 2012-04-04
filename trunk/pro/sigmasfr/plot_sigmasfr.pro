pro plot_sigmasfr
; jm12mar19ucsd - make the sigma-sfr vs mass plot
    path = getenv('HIZEA_DATA')+'/sigmasfr/'

    hst = rsex(path+'hst_sample.dat')
    kcorr = mrdfits(path+'sigmasfr_kcorrect.fits.gz',1)
    lir = mrdfits(path+'sigmasfr_lir.fits.gz',1)

    veilleux = rsex(path+'06veilleux.dat')

    ;xrange = [-22,-26]
    xrange = [8.3, 12.2]
    yrange = [-1.5,4.8]
    ;sigmasfr_edd = alog10(5000)
    sigmasfr_edd = alog10(2000)    
    sigmasfr_meurer = alog10(45.)    

    psfile = path+'sigmasfr.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=2, $
      xmargin=[1.5,0.4], width=6.6, height=5.3
    djs_plot, [0], [0], /nodata, position=pos, xsty=5, ysty=5, $
      xrange=xrange, yrange=yrange
    ;polyfill, [xrange,reverse(xrange)], [xrange*0+sigmasfr_edd,xrange*0+yrange[1]], $
    ;  /data, color=im_color('grey90'), noclip=0, /line_fill, orientation=135, spacing=0.15
    ;polyfill, [xrange,reverse(xrange)], [xrange*0+sigmasfr_edd,xrange*0+yrange[1]], $
    ;  /data, color=im_color('grey90'), noclip=0, /line_fill, orientation=45, spacing=0.15
    djs_oplot, !x.crange, sigmasfr_edd*[1,1], line=0, color=im_color('grey90')
    djs_oplot, !x.crange, sigmasfr_meurer*[1,1], line=2, color=im_color('grey90')
    djs_plot, [0], [0], /noerase, position=pos, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange, $
      xtickinterval=1, xtitle='log (M_{*} / M'+sunsymbol()+')', $;xtitle='Absolute H Magnitude', $
      ytitle='log (\Sigma_{SFR} / M'+sunsymbol()+' yr^{-1} kpc^{-2})'

    ;im_legend, 'Eddington-limited Star Formation', box=0, $
    ;  charsize=1.2, pos=[8., 3.7]
    im_legend, 'Eddington limit', box=0, $
      charsize=1.2, pos=[11.15, sigmasfr_edd+0.33]
    im_legend, 'Meurer+97 limit', box=0, $
      charsize=1.2, pos=[11.15, sigmasfr_meurer+0.33]

    ;im_legend, ['HST-WISE Sample','Veilleux+06','Law+12','Overzier+09'], $
    im_legend, [textoidl('HST-WISE Sample (z\sim0.6, symsize\propto|v_{out}|)'), $
    ;symsize\propto\r_{e}), symsize\propto|v_{out}|
      textoidl('gas-rich mergers (z<0.3)'), $
      textoidl('Lyman break galaxies (z=1.5-3.6)'), $
      textoidl('Lyman break analogs (z<0.3)')], $
      pos=[xrange[0]+0.1, yrange[1]-0.1], box=0, $
      charsize=1.5, psym=[16,15,14,46], color=['black','orange','blue','dark green']
    
; HST sample
    good = where(lir.sfr_chary gt -900)    
    sigmasfr = alog10(lir[good].sfr_chary/(!pi*hst[good].r_e^2)/2.0)
    absmag_h = kcorr[good].k_fnuvugrizjhk_absmag_00[8]
    mtest = (4.83-(absmag_h))/2.5
    ml = 10.^(kcorr[good].k_mass) / 10.^(mtest)
    masstolight = median(ml)
;   niceprint, absmag_h, sigmasfr
    ;djs_oplot, absmag_h, sigmasfr, psym=symcat(16), symsize=1.5

;; point size based on effective radius
;    maxre = max(alog10(hst.r_e))
;    minre = min(alog10(hst.r_e))
;    rearr = findgen(101) / (100/(maxre-minre)) + minre
;    maxsize = 2.0
;    minsize = 1.0
;    sizearr = findgen(101) / (100/(maxsize-minsize)) + minsize
;    for i=0L, n_elements(good)-1L do begin      
;      psize = spline(rearr, sizearr, alog10(hst[good[i]].r_e))
;      djs_oplot, [0., kcorr[good[i]].k_mass], [0., sigmasfr[i]], psym=symcat(16), symsize=psize    
;    endfor

; point size based on velocity
    outflow = where(abs(hst.vout) gt 0., noutflow)
    maxout = alog10(max(abs(hst[outflow].vout)))
    minout = alog10(min(abs(hst[outflow].vout)))
    outarr = findgen(101) / (100./(maxout-minout)) + minout
    maxsize = 2.0
    minsize = minout / maxout * maxsize
    ;minsize = 1.0
    minsize = 1.4
    sizearr = findgen(101) / (100/(maxsize-minsize)) + minsize
    for i=0L, n_elements(good)-1L do begin     
      if abs(hst[good[i]].vout) gt 0. then $ 
        psize = spline(outarr, sizearr, alog10(abs(hst[good[i]].vout))) else $
        psize = 1.0
      djs_oplot, [0., kcorr[good[i]].k_mass], [0., sigmasfr[i]], psym=symcat(16), symsize=psize    
    endfor

; Veilleux+06
    sigmasfr = alog10(4.5D-44*10.0^veilleux.lir*im_lsun()/(!pi*veilleux.r50^2)/2.0)
    ;v06_mass = (4.83-(veilleux.absmag_h))/2.5 ; assuming M/L = 1 
    v06_mass = (4.83-(veilleux.absmag_h))/2.5 + alog10(masstolight) ; assuming median M/L of HST sample 
    ;djs_oplot, veilleux.absmag_h, sigmasfr, psym=symcat(15), color='orange'
    djs_oplot, v06_mass, sigmasfr, psym=symcat(15), color='orange'

; Law+12
     law = rsex(path+'12law.dat')
     oplot, alog10(law.mstar), alog10(law.sigma_sfr / 2.), psym=symcat(14), color=djs_icolor('blue')     

; Overzier+09
     over = rsex(path+'09overzier.dat')
     oplot, over.mgal, alog10(over.sfr / (2.*!pi*over.re^2)), $
       psym=symcat(46), color=djs_icolor('dark green'), symsize=2.0    

    im_plotconfig, psfile=psfile, /psclose, /pdf, /pskeep


stop
return
end
