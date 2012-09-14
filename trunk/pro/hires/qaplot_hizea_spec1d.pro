pro qaplot_hizea_spec1d, datapath=datapath, qafile=qafile
; jm09dec10ucsd - build a QAplot for each individual 1D spectrum

    if (n_elements(datapath) eq 0) then datapath = './'
    im_plotconfig, 2, pos, psfile=qafile, yspace=0.6, $
      xspace=0.1, width=4.5*[1,1], height=3.2*[1,1], $
      charsize=1.6, ymargin=[0.5,1.0]

    fits = file_search(datapath+'spec1d/*.fits.gz',count=nfits)

; zoom into five spectral regions:
;  Mg II - 2796,2803
;  Mg I  - 2852
;  Fe II - 2344,2374,2382
;  Fe II - 2586,2600 *and* Mn II 2576,2594,2606

    linewave = [2796,2803,2852,2344,2374,2382,$
      2586,2600,2576,2594,2606]
    
    xrange = fltarr(2,4)
    xrange[*,0] = minmax([2344,2374,2382])+5*[-1,1]
    xrange[*,1] = minmax([2586,2600,2576,2594,2606])+5*[-1,1]
    xrange[*,2] = minmax([2796,2803])+[-25,+15]
    xrange[*,3] = 2852+[-25,+15]

    linelabel = [$
      'FeII \lambda\lambda2344,2374,2382',$
      'FeII \lambda\lambda2586,2600, MnII \lambda\lambda2576,2594,2606',$
      'MgII \lambda\lambda2796,2803',$
      'MgI \lambda2852']
    
    colors = ['orange','blue','red','green','purple','cyan','yellow']
    pad = 10
    for ii = 0, nfits-1 do begin
       splog, 'Reading '+fits[ii]
       spec = mrdfits(fits[ii],1,/silent)
       norder = n_elements(spec.order)

       zobj = spec.z
       wave = spec.wave/(1.0+zobj)
       srt = sort(wave)
       wave = wave[srt]
       flux = spec.flux[srt]
       ferr = spec.ferr[srt]
       snr = flux/(ferr+(ferr eq 0))*(ferr ne 0)

; zoom into the lines of interest
       for jj = 0, 3 do begin
          xrange1 = xrange[*,jj]
; set the y-ranges according to Mg II          
          get_element, wave, xrange1, xx
;         yrange1 = [0.0,im_max(flux[xx[0]:xx[1]],sigrej=3.0)*1.05]
;         yrange1 = [0.0,im_max(snr[xx[0]:xx[1]],sigrej=3.0)*1.2]
          yrange1 = [0,13]

          if (odd(jj) eq 0) then begin
;            ytitle1 = 'Flux (counts s^{-1})'
             ytitle1 = 'S/N (pixel^{-1})'
             delvarx, ytickname
          endif else begin
             ytitle1 = ''
             ytickname = replicate(' ',10)
          endelse
          
          if (jj le 1) then xtitle1 = '' else $
            xtitle1 = 'Rest Wavelength (\AA)'
          
          djs_plot, [0], [0], /nodata, xsty=1, ysty=3, xrange=xrange1, $
            yrange=yrange1, ytickname=ytickname, xtitle=xtitle1, $
            ytitle=ytitle1, position=pos[*,jj], noerase=(jj gt 0)
          for ll = 0, n_elements(linewave)-1 do djs_oplot, linewave[ll]*[1,1], $
            !y.crange, line=1
;         djs_oplot, !x.crange, [0,0], line=0, thick=2
          ccount = 0
          for kk = 0, norder-1 do begin
             good = where(spec.wave[*,kk] gt 0.0,ngood)
             orderwave = spec.wave[good[pad:ngood-1-pad],kk]/(1.0+zobj)
             orderflux = spec.flux[good[pad:ngood-1-pad],kk]
             orderferr = spec.ferr[good[pad:ngood-1-pad],kk]
             if (max(orderwave) gt xrange1[0]) and $
               (min(orderwave) lt xrange1[1]) then begin
;               djs_oplot, orderwave, orderflux/orderferr, color=colors[ccount]
                djs_oplot, orderwave, smooth(orderflux/orderferr,3), $
                  color=colors[ccount]
                ccount++
             endif
          endfor 
          im_legend, linelabel[jj], /left, /top, box=0, charsize=1.2, margin=0
       endfor
       thistitle = strtrim(spec.galaxy,2)+', exp'+$
         string(spec.frame,format='(I4.4)')+$
         ', z='+strtrim(string(spec.z,format='(F12.5)'),2)
       xyouts, pos[2,0], pos[3,0]+0.015, thistitle, align=0.5, /norm, $
         charsize=1.4
    endfor
       
;       djs_plot, [0], [0], /nodata, /noerase, xsty=1, ysty=1, $
;         xrange=xrange2, yrange=yrange2, xtitle='Rest Wavelength (\AA)', $
;         ytitle='S/N (pixel^{-1})', position=pos[*,1]
;       djs_oplot, !x.crange, [0,0], line=0, thick=2
;       ccount = 0
;       for jj = 0, norder-1 do begin
;          good = where(spec.wave[*,jj] gt 0.0,ngood)
;          wave = spec.wave[good[pad:ngood-1-pad],jj]/(1.0+zobj)
;          flux = spec.flux[good[pad:ngood-1-pad],jj]
;          ferr = spec.ferr[good[pad:ngood-1-pad],jj]
;          if (max(wave) gt xrange2[0]) and (min(wave) lt xrange2[1]) then begin
;             djs_oplot, wave, smooth(flux/ferr,4), color=colors[ccount]
;             ccount++
;          endif
;       endfor
;    endfor
    im_plotconfig, psfile=qafile, /psclose, /gzip

return
end
