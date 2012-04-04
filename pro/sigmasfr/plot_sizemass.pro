pro plot_sizemass
;+
;  amd120330 -- make a size-mass plot for sigmasfr paper
;
;
;
;-

    path = getenv('HIZEA_DATA')+'/sigmasfr/'
    hst = rsex(path+'hst_sample.dat')
    kcorr = mrdfits(path+'sigmasfr_kcorrect.fits.gz',1)
    lir = mrdfits(path+'sigmasfr_lir.fits.gz',1)

    xrange = [10.3, 11.7]
    yrange = [-1.2, 1.5]
    xtitle = 'log (M_{*} / M'+sunsymbol()+')'
    ytitle=textoidl('log (r_{e} / kpc)')

    psfile = path+'sizemass_sigmasfr.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=2, $
      xmargin=[1.5,0.4], width=6.6, height=5.3
    djs_plot, [0], [0], xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle

; hst sample
    good = where(lir.sfr_chary gt -900)    
    sigmasfr = alog10(lir[good].sfr_chary/(!pi*hst[good].r_e^2)/2.0)
    absmag_h = kcorr.k_fnuvugrizjhk_absmag_00[8]
    mtest = (4.83-(absmag_h))/2.5
    ml = 10.^(kcorr.k_mass) / 10.^(mtest)
    masstolight = median(ml)
; point size based on sigmasfr
    maxsig = max(sigmasfr)
    minsig = min(sigmasfr)
    sigarr = findgen(101) / (100/(maxsig-minsig)) + minsig
    maxsize = 2.0
    minsize = 1.0
    sizearr = findgen(101) / (100/(maxsize-minsize)) + minsize
    for i=0L, n_elements(good)-1L do begin      
      psize = spline(sigarr, sizearr, sigmasfr[i])
      djs_oplot, [0., kcorr[good[i]].k_mass], [0., alog10(hst[good[i]].r_e)], psym=symcat(16), symsize=psize    
    endfor
    djs_oplot, kcorr.k_mass, alog10(hst.r_e), psym=symcat(16), symsize=1.0
;; point size based on velocity
;    outflow = where(abs(hst.vout) gt 0., noutflow)
;    maxout = alog10(max(abs(hst[outflow].vout)))
;    minout = alog10(min(abs(hst[outflow].vout)))
;    outarr = findgen(101) / (100./(maxout-minout)) + minout
;    maxsize = 2.0
;    minsize = 1.4
;    sizearr = findgen(101) / (100/(maxsize-minsize)) + minsize
;    for i=0L, n_elements(hst)-1L do begin     
;      if abs(hst[i].vout) gt 0. then $ 
;        psize = spline(outarr, sizearr, alog10(abs(hst[i].vout))) else $
;        psize = 1.0
;      djs_oplot, [0., kcorr[i].k_mass], [0., alog10(hst[i].r_e)], psym=symcat(16), symsize=psize    
;    endfor


  
; guo et al. 2009    
    maxis = im_array(10.0,12.5,0.02)
    djs_oplot, maxis, poly(maxis,[-8.45,0.83]), line=5, thick=5


; hopkins et al. 2010 sigma_mass
    djs_oplot, maxis, alog10(sqrt(10^(maxis)/(2*!pi*1.d11))), linestyle=1
    ;djs_oplot, maxis, alog10(sqrt(10^(maxis)/(2*!pi*1.d12))), linestyle=1


;; veilleux et al. 2006    
    veilleux = rsex(path+'06veilleux.dat')
  ; assuming median M/L of HST sample 
    v06_mass = (4.83-(veilleux.absmag_h))/2.5 + alog10(masstolight) 
    djs_oplot, v06_mass, alog10(veilleux.r50), psym=symcat(15), color='orange'

; van dokkum et al. 2008
    v08_re = [0.47,0.49,0.78,0.76,0.92,0.93,1.42,1.89,2.38]
    v08_mass = 1D11*[0.6,0.9,2.5,3.0,1.4,2.0,1.8,1.7,2.0]

    psize = 1.5

    plotsym, 3, psize, /fill, color=djs_icolor('red') 
    djs_oplot, alog10(v08_mass), alog10(v08_re), psym=8

;; Overzier+09 ; issue of clump mass v. galaxy mass
;     over = rsex(path+'09overzier.dat')
;     ;oplot, over.mgal, alog10(over.re), $
;     oplot, over.mclump, alog10(over.re), $
;       psym=symcat(46), color=djs_icolor('dark green'), symsize=2.0  

; trujillo et al. 2007
  xmass = 11.2
  z0size = 0.84
  ;oplot, [xmass, 0.], [z0size, 0.], psym=6
  xoff = 0.025
  yoff = -0.03
  tcolor = 'brown'
  tthick = 8.
  zsize=1.8
  ;zstr = ['0.1', '0.2<z<0.5', '0.5<z<0.8', '0.8<z<1.1', '1.1<z<1.4', $
  ;        '1.4<z<1.7', '1.7<z<2.0', '2.0<z<2.5', '2.5<z<3.0']
  ;fac = [1., 0.84, 0.63, 0.41, 0.34, 0.26, 0.23, 0.23, 0.14]
  xyouts, xmass+xoff, z0size+alog10(0.24)+yoff, textoidl('z\sim2.0'), $
    charsize=zsize, align=0.5, color=djs_icolor(tcolor), charthick=tthick
  xyouts, xmass+xoff, z0size+alog10(0.40)+yoff, textoidl('z\sim1.0'), $
    charsize=zsize, align=0.5, color=djs_icolor(tcolor), charthick=tthick
  xyouts, xmass+xoff, z0size+alog10(0.73)+yoff, textoidl('z\sim0.5'), $
    charsize=zsize, align=0.5, color=djs_icolor(tcolor), charthick=tthick
  ;xyouts, xmass+xoff, z0size+yoff, textoidl('z\sim0.0'), $
  ;  charsize=zsize, align=0.5, color=djs_icolor(tcolor), charthick=tthick
  t07_mass = [xmass, xmass, xmass]
  t07_size = z0size+alog10([0.24, 0.40, 0.73])
  ;djs_oplot, t07_mass, t07_size, psym=symcat(14), symsize=psize, $
  ;  color=tcolor, thick=10.
  plots, [xmass, xmass]+0.15, [t07_size[0], z0size], linestyle=0, $
    thick=tthick*2, color=djs_icolor(tcolor)
  djs_oplot, [0., xmass+0.15], [0., z0size], psym=symcat(17), $
    color=djs_icolor(tcolor), symsize=psize*2, thick=tthick*2

    xyouts, xmass+0.1, z0size+0.12, $
      textoidl('Size-Mass Relation at z\sim0.1'), $
      orientation=19, align=0.5, charsize=1.3
    xyouts, xmass+0.2, t07_size[0]+0.03, 'Trujillo+07', $;'Size Evolution', $
      orientation=90, color=djs_icolor(tcolor), charsize=1.3, $
      charthick=6.
    xyouts, xmass+0.25, z0size-1.14, $
      textoidl('\Sigma_{*}=10^{11} M')+sunsymbol()+textoidl(' kpc^{-2}'), $
      orientation=12, align=0.5, charsize=1.5
    ;xyouts, xmass+0.25, z0size-1.64, $
    ;  textoidl('\Sigma_{*}=10^{12} M')+sunsymbol()+textoidl(' kpc^{-2}'), $
    ;  orientation=12, align=0.5, charsize=1.5

  lpos = [xrange[0]+0.05, yrange[1]-0.05]

  im_legend, [textoidl('HST-WISE Sample (z\sim0.6, symsize\propto\Sigma_{SFR})'), $
; symsize\propto\Sigma_{SFR}, symsize\propto|v_{out}|
      textoidl('compact quiescent galaxies (z\sim2.3)'), $
      textoidl('gas-rich mergers (z<0.3)')], $
      pos=lpos, box=0, $
      charsize=1.5, psym=[16,3,15], $
      color=['black','red','orange']
   
   plotsym, 3, psize*1.2, /fill, color=djs_icolor('red') 
   ;djs_oplot, [0., 10.453], [0., 0.865], psym=8
   djs_oplot, [0., lpos[0]+0.053], [0., lpos[1]-0.226], psym=8
   

    im_plotconfig, psfile=psfile, /psclose, /pdf, /pskeep


stop

end
