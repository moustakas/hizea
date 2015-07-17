pro sizemass

; van dokkum et al. 2010    
    mass = [11.45,11.36,11.28,11.21,11.15]
    re = [12.4,8.0,5.3,4.1,3.0]
    re_err = [1.6,1.2,0.3,0.3,0.4]

; van dokkum et al. 2008
    v08_re = [0.76,0.92,0.78,1.89,0.93,1.42,0.47,2.38,0.49]
    v08_re_err = [0.06,0.18,0.17,0.15,0.04,0.35,0.03,0.11,0.02]

    v08_re = [0.47,0.49,0.78,0.76,0.92,0.93,1.42,1.89,2.38]
    v08_re_err = [0.03,0.02,0.17,0.06,0.18,0.04,0.35,0.15,0.11]
    v08_mass = 1D11*[0.6,0.9,2.5,3.0,1.4,2.0,1.8,1.7,2.0]
    v08_mass_err = v08_mass*0.0+0.05
    
; guo et al. 2009    
    im_plotconfig, 0, pos, psfile='sizemass.eps', $
      xmargin=[1.4,0.3], width=6.8, height=6.5
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xrange=[10.5,11.6], yrange=[-0.75,1.23], $
      xtitle='log_{10} (M/M'+sunsymbol()+')', ytitle='log_{10} (r_{e}) (h_{70}^{-1} kpc)'
;; van dokkum et al. 2010    
;    oploterror, mass, alog10(re), re_err/re/alog(10.0), $
;      psym=symcat(6,thick=8), symsize=2.0, color=fsc_color('dodger blue',101), $
;      errcolor=fsc_color('dodger blue',101), errthick=8
;    label = ['z=0','z=0.6','z=1.1','z=1.6','z=2']
;    align = [1.0,0.0,0.0,0.0,0.0]
;    moff = [-0.05,0.05,0.05,0.05,0.05]
;    for ii = 0, n_elements(mass)-1 do xyouts, mass[ii]+moff[ii], $
;      alog10(re[ii])-0.03, label[ii], /data, align=align[ii], $
;      charsize=1.5, color=fsc_color('dodger blue',101)
;;   xyouts, 11.35, 0.35, 'van Dokkum+10', align=0.5, charsize=1.3
; van dokkum et al. 2008
    oploterror, alog10(v08_mass), alog10(v08_re), v08_mass_err, $
      v08_re_err/v08_re/alog(10.0), psym=symcat(9,thick=8), symsize=1.6, $
      color=fsc_color('orange',103), $
      errcolor=fsc_color('orange',103), errthick=8
    xyouts, 11.3, 0.5, 'van Dokkum et al. 2008 (z~2.3)', align=0.5, $
      charsize=1.3, color=fsc_color('orange',103)
;   plots, alog10(1.7D11), alog10(0.9), psym=symcat(9,thick=12), $
;     symsize=4.0, color=fsc_color('orange',103)
; guo et al. 2009    
    maxis = im_array(10.0,12.5,0.02)
    djs_oplot, maxis, poly(maxis,[-8.45,0.83]), line=5, thick=5
    xyouts, 11.05, 0.77, 'Size-Mass Relation at z~0.1', $
      orientation=24, align=0.5, charsize=1.5
; our sample
    polyfill, [10.57,10.57,11.27,11.27,10.57], $
      [-0.7,0.1,0.1,-0.7,-0.7], /data, /line_fill, noclip=0, $
      color=fsc_color('tan',102), orientation=45, spacing=0.1
    xyouts, 10.92, -0.2, 'Post-Starbursts at z=0.4-0.8', $
      align=0.5, charsize=1.5
    xyouts, 11.3, -0.32, '?', align=0.5, charsize=2.5
    arrow, 11.3, -0.2, 11.3, 0.1, /data, hsize=-0.15, hthick=8, thick=8
    arrow, 11.3, -0.37, 11.3, -0.7, /data, hsize=-0.15, hthick=8, thick=8
; recent HST observations
    plots, 10.69, alog10(6.54*1.16*0.04), psym=symcat(4,thick=12), symsize=4.0, $
      color=fsc_color('firebrick',100)
    plots, 10.72, alog10(6.70*1.09*0.04), psym=symcat(4,thick=12), symsize=4.0, $
      color=fsc_color('firebrick',100)
    xyouts, 10.75, -0.49, 'J1104+5946 (z=0.57)', /data, align=0.0, charsize=1.5
    xyouts, 10.76, -0.57, 'J0826+4305 (z=0.60)', /data, align=0.0, charsize=1.5
; legend
    label = ['Quiescent Galaxies at z~2.3','Average Galaxies at z=0-2']
    psym = [9,6]
    color = ['orange','dodger blue']
;   im_legend, label, /left, /top, box=0, charsize=1.3, $
;     psym=psym, symsize=2.0, color=color, thick=12.0, $
;     charthick=2.5
    im_plotconfig, /psclose
    
return
end
    
