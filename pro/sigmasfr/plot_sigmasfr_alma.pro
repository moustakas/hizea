pro plot_sigmasfr_alma
; jm12mar19ucsd - make the sigma-sfr vs mass plot
    path = getenv('HIZEA_DATA')+'/sigmasfr/'

    hst = rsex(path+'hst_sample.dat')
    kcorr = mrdfits(path+'sigmasfr_kcorrect.fits.gz',1)
    lir = mrdfits(path+'sigmasfr_lir.fits.gz',1)
    isedfit = mrdfits(path+'isedfit/sigmasfr_fsps_chab_smc_sfhgrid01.fits.gz', 1)

    veilleux = rsex(path+'06veilleux.dat')

    ;xrange = [-22,-26]
    xrange = [8.3, 12.0]
    ;xrange = [8.6, 12.0]
    ;xrange = [9.0, 12.0]
    ;xrange = [9.9, 11.6]
    ;yrange = [-1.5,4.8]
    yrange = [-2.5,5.5]
    ;sigmasfr_edd = alog10(5000)
    ;sigmasfr_edd = alog10(2000)
    sigmasfr_edd = alog10(3000.)    
    sigmasfr_meurer = alog10(45.)-0.25    
    sigmasfr_heckman = alog10(0.1)    

  ;!p.thick=3
  ;!x.thick=3
  ;!y.thick=3
  ;!p.charthick=3

    psfile = path+'sigmasfr_alma.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=2, $
      xmargin=[1.5,0.4], width=6.6, height=5.3
    djs_plot, [0], [0], /nodata, position=pos, xsty=5, ysty=5, $
      xrange=xrange, yrange=yrange
    ;polyfill, [xrange,reverse(xrange)], [xrange*0+sigmasfr_edd,xrange*0+yrange[1]], $
    ;  /data, color=im_color('grey90'), noclip=0, /line_fill, orientation=135, spacing=0.15
    ;polyfill, [xrange,reverse(xrange)], [xrange*0+sigmasfr_edd,xrange*0+yrange[1]], $
    ;  /data, color=im_color('grey90'), noclip=0, /line_fill, orientation=45, spacing=0.15
    ;djs_oplot, !x.crange, sigmasfr_edd*[1,1], line=0, color=im_color('grey90')
    djs_oplot, [11.0, 12.0], sigmasfr_edd*[1,1], line=0, color=djs_icolor('black')
    ;djs_oplot, !x.crange, sigmasfr_meurer*[1,1], line=2, color=im_color('grey90')
    djs_oplot, [11.0, 12.0], sigmasfr_meurer*[1,1], line=2
    ;djs_oplot, !x.crange, sigmasfr_heckman*[1,1], line=1, color=im_color('grey90')
    djs_oplot, [11.0, 12.0], sigmasfr_heckman*[1,1], line=1

    djs_plot, [0], [0], /noerase, position=pos, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange, $
      xtickinterval=1, xtitle='log (M_{*} / M'+sunsymbol()+')', $;xtitle='Absolute H Magnitude', $
      ytitle='log (\Sigma_{SFR} / M'+sunsymbol()+' yr^{-1} kpc^{-2})'

    ;im_legend, 'Eddington-limited Star Formation', box=0, $
    ;  charsize=1.2, pos=[8., 3.7]
    yoff = 0.45
    im_legend, 'Eddington limit', box=0, $
      charsize=1.2, pos=[11.08, sigmasfr_edd+yoff]
    im_legend, 'Meurer+97 limit', box=0, $
      charsize=1.2, pos=[11.03, sigmasfr_meurer+yoff]
    im_legend, 'threshold for winds', box=0, $
      charsize=1.2, pos=[10.94, sigmasfr_heckman+yoff]

    ;lpos = [xrange[0]+0.1, yrange[1]-0.1]
    ;lpos = [xrange[0]+0.5, yrange[1]-0.2]

;    ;im_legend, ['HST-WISE Sample','Veilleux+06','Law+12','Overzier+09'], $
;    im_legend, [textoidl('HST-WISE Sample (z\sim0.6, symbol size\proptov_{out})'), $
;    ;symsize\propr_{e}), symsize\propto|v_{out}|
;      textoidl('gas-rich mergers (z<0.3)')] , $
;      ;textoidl('Lyman break galaxies (z=1.5-3.6)'), $
;      ;textoidl('Arp 220 nucleus'), $
;      ;textoidl('Lyman break analogs (z<0.3)'), $
;      ;textoidl('star-forming galaxies (z=1.5-3.6)')], $
;      pos=lpos, box=0, $
;      spacing=2.5, $
;      ;charsize=1.4, psym=[16,15,17,46,14], color=['black','orange','magenta', 'dark green','blue']
;      charsize=1.4, psym=[16,15], color=['black','orange']

    ;lpos = [10.15, 5.25]
    lpos = [9.55, 5.25]
    im_legend, textoidl('HST-WISE sample (this paper, z\sim0.6)'), $;, symbol size\proptov_{out})'), $
      pos=lpos, box=0, charsize=1.4, psym=16
    ;arrow, 10.7, 5.0, 11.0, 4.0, /data, thick=5, hsize=-0.2
    arrow, 11.0, 4.8, 11.0, 4.0, /data, thick=5, hsize=-0.3
    
    mpos = [8.95, 4.45]
    ;mpos = [9.0, 4.45]
    im_legend, textoidl('gas-rich mergers (z<0.3)'), $
      pos=mpos, box=0, charsize=1.4, psym=15, color='orange'
    ;arrow, 10.0, 3.9, 10.3, 2.9, /data, thick=5, hsize=-0.2
    arrow, 10.1, 3.9, 10.3, 3.1, /data, thick=5, hsize=-0.28

    
    ;xyouts, xrange[0]+0.3, yrange[0]+0.5, 'star-forming galaxies, 50% contours: ', charsize=1.5
    ;xyouts, xrange[0]+0.3, yrange[0]+0.3, 'star-forming galaxies, 68% contours: ', charsize=1.5
    ;xyouts, xrange[0]+0.3, yrange[1]-1.3, 'star-forming galaxies, 68% contours: ', charsize=1.4
    ;xyouts, xrange[0]+0.1, yrange[1]-2.1, 'star-forming galaxies', charsize=1.4
    ;xyouts, xrange[0]+0.1, yrange[1]-2.3, '68% contours:', charsize=1.4

    spos = [8.3, 3.4]
    ;spos = [8.3, -0.5]
    im_legend, textoidl('star-forming galaxies (z\sim1)'), $
      pos=spos, box=0, charsize=1.4 
    ;arrow, 9.1, 3.1, 9.4, 2.1, /data, thick=5, hsize=-0.2
    arrow, 9.2, 2.8, 9.2, 2.0, /data, thick=5, hsize=-0.3


    ;zpos = [xrange[1]-1.1, yrange[0]+1.7]
    zpos = [xrange[0], yrange[1]-2.2]
    

    im_legend, '68%', pos=[9., 0.], box=0, charsize=1.4, color=djs_icolor('red')
    im_legend, '95%', pos=[9., 0.8], box=0, charsize=1.4, color=djs_icolor('red')
    im_legend, '99.7%', pos=[9., 1.9], box=0, charsize=1.4, color=djs_icolor('red')

    ;im_legend, ['68%', '95%', '99.7%'], pos=zpos, charsize=1.5, box=0, $
    ;  linestyle=[0,2,1], color=[djs_icolor('red'), djs_icolor('red'), djs_icolor('red')]

   ; im_legend, [textoidl('z\sim3'), textoidl('z\sim2'), $
   ;   textoidl('z\sim1'), $
   ;   textoidl('z\sim0.1')], linestyle=[3, 1, 2, 0], $
   ;   color=[djs_icolor('blue'), djs_icolor('purple'), djs_icolor('red'), djs_icolor('black')], $
   ;   ;/bottom, /right, $
   ;   pos=zpos, $
   ;   charsize=1.5, box=0, thick=[10., 10., 10., 5.]

; HST sample
    good = where(lir.sfr_chary gt -900)    
    sigmasfr = alog10(lir[good].sfr_chary/(!pi*hst[good].r_e^2)/2.0)
  ; convert to chabrier IMF
    sigmasfr = sigmasfr - 0.25
    absmag_h = kcorr[good].k_fnuvugrizjhk_absmag_00[8]
    mtest = (4.83-(absmag_h))/2.5
    ml = 10.^(kcorr[good].k_mass) / 10.^(mtest)
    masstolight = median(ml)
;   niceprint, absmag_h, sigmasfr
    ;djs_oplot, absmag_h, sigmasfr, psym=symcat(16), symsize=1.5
    ;djs_oploterr, [0., 11.5], [0., 2.5], yerr=[0., 0.301], /cap
    oploterror, [0., 11.5], [0., 4.3], [0., 0.2], [0., 0.3], psym=3

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



; Veilleux+06
    v06_sigmasfr = alog10(4.5D-44*10.0^veilleux.lir*im_lsun()/(!pi*veilleux.r50^2)/2.0)-0.25
    ;v06_mass = (4.83-(veilleux.absmag_h))/2.5 ; assuming M/L = 1 
    v06_mass = (4.83-(veilleux.absmag_h))/2.5 + alog10(masstolight) ; assuming median M/L of HST sample 
    ;djs_oplot, veilleux.absmag_h, sigmasfr, psym=symcat(15), color='orange'
    djs_oplot, v06_mass, v06_sigmasfr, $
      psym=symcat(15), color='orange'
      ;psym=symcat(46), color=djs_icolor('dark green'), symsize=2.0 
      ;psym=symcat(14), color=djs_icolor('blue'), symsize=1.8
; this paper
    for i=0L, n_elements(good)-1L do begin     
      if abs(hst[good[i]].vout) gt 0. then $ 
        psize = spline(outarr, sizearr, alog10(abs(hst[good[i]].vout))) else $
        psize = 1.0
      ;djs_oplot, [0., kcorr[good[i]].k_mass], [0., sigmasfr[i]], $
      ;  psym=symcat(16), symsize=psize    
      djs_oplot, [0., isedfit[good[i]].mass], [0., sigmasfr[i]], $
        ;psym=symcat(16), symsize=psize    
        psym=symcat(9), symsize=psize
      if hst[good[i]].z lt 0.64 and hst[good[i]].dec lt 35. then $
      djs_oplot, [0., isedfit[good[i]].mass], [0., sigmasfr[i]], $
        psym=symcat(16), symsize=psize    
        ;psym=symcat(9), symsize=psize      

    endfor

;; Law+12
     law = rsex(path+'12law_v2.dat')
     ;oplot, alog10(law.mstar), alog10(law.sigma_sfr / 2.), psym=symcat(14), color=djs_icolor('blue'), $
     ;  symsize=psize/2  
     ;highz = where(law.zspec gt 2.5, nhighz)   
     ;oplot, alog10(law[highz].mstar), alog10(law[highz].sigma_sfr / 2.), psym=symcat(14), color=djs_icolor('blue'), $
     ;  symsize=psize


; Overzier+09
     over = rsex(path+'09overzier.dat') ; mgal v. mclump
     ; convert kroupa to chabrier imf
     djs_oplot, over.mgal, alog10(over.sfr / (2.*!pi*over.re^2)) - 0.074, $
       psym=symcat(15), color='orange'
       ;psym=symcat(46), color=djs_icolor('dark green'), symsize=2.0    

; Arp 220, Sigma_SFR from Ken98 (2.98) 
; stellar mass from Rodriguez-Zaurin+08 is 3.d10
; stellar mass from Scoville+97 is 2.5d9 (just for the nucleus)
    ;oplot, [0., alog10(2.5d9)], [0., 2.98], psym=symcat(17), color=djs_icolor('magenta'), symsize=psize
    ;xyouts, alog10(2.5d9)+0.06, 2.98-0.06, 'Arp 220', charsize=1.2
    djs_oplot, [0., alog10(3.0d10)], [0., 2.98-0.25], $
      ;psym=symcat(17), color=djs_icolor('magenta'), symsize=psize
      psym=symcat(15), color='orange'
   ;xyouts, alog10(3.0d10)-0.45, 2.98-0.08, 'Arp 220', charsize=1.2

; Rubin+10
   rubin = rsex(path+'10rubin.dat')
   ;djs_oplot, rubin.mstar, alog10(10^(rubin.sfr)/(2.*!pi*rubin.re^2)), psym=symcat(46), color='dark green'

; Wuyts+11
    z0 = rsex(path+'wuyts/logM_logSFRsurfdens_0.02z0.20.cat')
    goodz0 = where(10^(z0.sfr) / 10.^(z0.mstar) gt 1. / 12.6d9, ngoodz0)
    z1 =  rsex(path+'wuyts/logM_logSFRsurfdens_0.5z1.5.cat')
    goodz1 = where(10^(z1.sfr) / 10.^(z1.mstar) gt 1. / 5.75d9, ngoodz1)
    z2 =  rsex(path+'wuyts/logM_logSFRsurfdens_1.5z2.5.cat')
    goodz2 = where(10^(z2.sfr) / 10.^(z2.mstar) gt 1. / 3.23d9, ngoodz2)



  ;oplot, z0.mstar, z0.sigmasfr, psym=3
  ;oplot, z1.mstar, z1.sigmasfr, psym=3, color=djs_icolor('blue')
  ;oplot, z2.mstar, z2.sigmasfr, psym=3, color=djs_icolor('red')

  ;levels = [0.68, 0.95]
  ;levels = [0.25, 0.75]
  ;levels = [0.75]
  ;levels = [0.50]
  ;levels = [0.90]
  levels = 0.68


  ;im_hogg_scatterplot, alog10(law[highz].mstar), alog10(law[highz].sigma_sfr / 2.), $
  ;  /noerase, xsty=5, xrange=xrange, $
  ;  yrange=yrange, ysty=5, pos=pos, /nogrey, /internal, $
  ;  xnpix=5, ynpix=5, levels=levels, $
  ;  contour_color=djs_icolor('blue'), $;('magenta blue'), $
  ;  c_linestyle=3, cthick=10

  ;im_hogg_scatterplot, z0[goodz0].mstar, z0[goodz0].sigmasfr, $
  ;  /noerase, xsty=5, xrange=xrange, $
  ;  yrange=yrange, ysty=5, pos=pos, /nogrey, /internal, $
  ;  xnpix=25., ynpix=25, $
  ;  ;levels=levels; cthick=2.
  ;  ;levels = [0.68, 0.95]
  ;  ;levels=[0.6,0.75,0.9,0.95]
  ;  ;levels=[0.25, 0.75]
  ;  levels=[0.68, 0.95, 0.997], /outliers

  im_hogg_scatterplot, z1[goodz1].mstar, z1[goodz1].sigmasfr, $
    /noerase, xsty=5, xrange=xrange, $
    yrange=yrange, ysty=5, pos=pos, /nogrey, /internal, $
    ;xnpix=25., ynpix=25, $
    xnpix=15, ynpix=15, $
    ;levels=levels, $
    levels=[0.68, 0.95, 0.997], $
    contour_color=djs_icolor('red'), $;('magenta blue'), $
    c_linestyle=[0,2,1], cthick=5;, /outliers

  ;im_hogg_scatterplot, z2[goodz2].mstar, z2[goodz2].sigmasfr, $
  ;  /noerase, xsty=5, xrange=xrange, $
  ;  yrange=yrange, ysty=5, pos=pos, /nogrey, /internal, $
  ;  xnpix=21., ynpix=21, $
  ;  ;levels=levels, $
  ;  levels=[0.68, 0.95, 0.997], $
  ;  contour_color=djs_icolor('purple'), $
  ;  c_linestyle=1, cthick=10, /outliers

    im_plotconfig, psfile=psfile, /psclose, /pdf, /pskeep


stop
return
end
