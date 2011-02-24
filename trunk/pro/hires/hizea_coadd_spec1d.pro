pro hizea_coadd_spec1d, hires1, info1, chip=chip, datapath=datapath, $
  qafile=qafile, fluxed=fluxed, std=std
; jm09sepucsd - coadd the 1D spectra, both individual orders and
;   repeat observations of the same object; if /FLUXED then coadd the
;   fluxed 1D spectra, otherwise coadd the counts

    if (n_elements(info1) ne n_elements(hires1)) then begin
       splog, 'Dimensions of INFO and HIRES do not match'
       return
    endif

    if (n_elements(datapath) eq 0L) then datapath = './'

; remove spectra that were never reduced (e.g., the red chip)
    for ii = 0, n_elements(info1)-1 do begin
       objfil = strtrim(hires1[ii].obj_fil,2)+'.gz'
       if file_test(datapath+objfil,/regular) then $
         if (n_elements(keep) eq 0) then keep = ii else $
           keep = [keep,ii]
    endfor
    info = info1[keep]
    hires = hires1[keep]

; condense the extracted spectra to a more convenient format;
; specifically, combine multiple chips of the same exposure/object
; into a single structure by grouping by unique galaxy+frame name
    allname = strtrim(info.galaxy,2)+'_'+string(info.frame,format='(I4.4)')
    uindx = uniq(allname,sort(allname))
    uname = allname[uindx]
    nobj = n_elements(uname)

    npix = 2010
    minnpix = 1900

    for ii = 0L, nobj-1L do begin
       these = where(uname[ii] eq allname,nthese)
       for kk = 0, nthese-1 do begin
          objfil = strtrim(hires[these[kk]].obj_fil,2)+'.gz'
          objstr = xmrdfits(datapath+objfil,1,STRUCTYP='hiresobjstrct',/silent)
          norder = n_elements(objstr)
          splog, 'Reading '+objfil+' norder = '+string(norder,format='(I0)')
          out = {chip: info[these[kk]].chip, order: -1, partial: 0, $
            npix: 0, wave: fltarr(npix), flux: fltarr(npix), ferr: fltarr(npix)}
          out = replicate(out,norder)
          out.order = objstr.order
          for jj = 0, norder-1 do begin
             good = where((objstr[jj].var gt 0.0),ngood)
             if (ngood lt minnpix) then begin
                out[jj].partial = 1 ; partial order
             endif else begin
                out[jj].npix = ngood
                out[jj].wave[0:ngood-1] = objstr[jj].wave[good]
                out[jj].flux[0:ngood-1] = objstr[jj].fx[good]/info[these[kk]].exp ; [counts/s]
                out[jj].ferr[0:ngood-1] = sqrt(objstr[jj].var[good])/info[these[kk]].exp
             endelse
          endfor
          if (kk eq 0) then allout = out else allout = [allout,out]
       endfor
; now recollapse ALLOUT into a one-element structure so that it can be
; combined with INFO; remove partial orders
       keep = where(allout.partial eq 0,norder)
       allout = allout[keep]
       final = {chip: intarr(norder), order: intarr(norder), npix: intarr(norder), $
         wave: fltarr(npix,norder), flux: fltarr(npix,norder), ferr: fltarr(npix,norder)}
       final.chip = allout.chip
       final.order = allout.order
       final.npix = allout.npix
       final.wave = allout.wave
       final.flux = allout.flux
       final.ferr = allout.ferr
       final = struct_addtags(struct_trimtags(info[these[0]],except='chip'),final)
; write out
       outfile = datapath+'spec1d/'+strtrim(info[these[0]].galaxy,2)+$
         '_'+string(info[these[0]].frame,format='(I4.4)')+'.fits'
       final.spec1dfile = file_basename(outfile)+'.gz'
       im_mwrfits, final, outfile
    endfor
    
;; make a QAplot
;    if (n_elements(qafile) ne 0L) then $
;      im_plotconfig, 9, pos, psfile=qafile, xmargin=[1.3,0.2]
;
;    fits = file_search(datapath+'spec1d/*.fits.gz',count=nfits)
;    
;    xrange1 = [3200,5300]       ; default range
;    xrange2 = [2750,2850]
;    colors = ['orange','blue','red','green','purple','cyan','yellow']
;    pad = 10
;    for ii = 0, nfits-1 do begin
;       splog, 'Reading '+fits[ii]
;       spec = mrdfits(fits[ii],1,/silent)
;       norder = n_elements(spec.order)
;; full observed wavelength range
;       good = where(spec.ferr gt 0.0)
;       flux = (spec.flux)[good]
;       ferr = (spec.ferr)[good]
;       snr = flux/ferr
;       yrange1 = [-0.005,im_max(flux,sigrej=3.0)*1.3]
;       yrange2 = [-0.005,im_max(snr,sigrej=3.0)*1.2]
;
;       djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=xrange1, $
;         yrange=yrange1, xtickname=replicate(' ',10), xtitle='', $
;         ytitle='Flux (counts s^{-1})', position=pos[*,0]
;       djs_oplot, !x.crange, [0,0], line=0, thick=2
;       for jj = 0, norder-1 do begin
;          good = where(spec.wave[*,jj] gt 0.0,ngood)
;          wave = spec.wave[good[pad:ngood-1-pad],jj]
;          flux = spec.flux[good[pad:ngood-1-pad],jj]
;          djs_oplot, wave, flux, psym=3
;       endfor
;       im_legend, [strtrim(spec.galaxy,2),'exp'+string(spec.frame,format='(I4.4)'),$
;         't_{exp}='+string(spec.exp,format='(I0)')+' s'], $
;         /left, /top, box=0, charsize=1.7
;
;       djs_plot, [0], [0], /nodata, /noerase, xsty=1, ysty=1, $
;         xrange=xrange1, yrange=yrange2, xtitle='Observed Wavelength (\AA)', $
;         ytitle='S/N', position=pos[*,1]
;       djs_oplot, !x.crange, [0,0], line=0, thick=2
;       for jj = 0, norder-1 do begin
;          good = where(spec.wave[*,jj] gt 0.0,ngood)
;          wave = spec.wave[good[pad:ngood-1-pad],jj]
;          flux = spec.flux[good[pad:ngood-1-pad],jj]
;          ferr = spec.ferr[good[pad:ngood-1-pad],jj]
;          djs_oplot, wave, flux/ferr, psym=3
;       endfor
;
;; zoom into Mg II
;       zobj = spec.z
;       wave = spec.wave/(1.0+zobj)
;       srt = sort(wave)
;       wave = wave[srt]
;       flux = spec.flux[srt]
;       ferr = spec.ferr[srt]
;       snr = flux/ferr
;       get_element, wave, xrange2, xx
;       yrange1 = [-0.009,im_max(flux[xx[0]:xx[1]],sigrej=3.0)*1.3]
;       yrange2 = [-0.009,im_max(snr[xx[0]:xx[1]],sigrej=3.0)*1.2]
;
;       djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=xrange2, $
;         yrange=yrange1, xtickname=replicate(' ',10), xtitle='', $
;         ytitle='Flux (counts s^{-1})', position=pos[*,0]
;       djs_oplot, !x.crange, [0,0], line=0, thick=2
;       ccount = 0
;       for jj = 0, norder-1 do begin
;          good = where(spec.wave[*,jj] gt 0.0,ngood)
;          wave = spec.wave[good[pad:ngood-1-pad],jj]/(1.0+zobj)
;          flux = spec.flux[good[pad:ngood-1-pad],jj]
;          if (max(wave) gt xrange2[0]) and (min(wave) lt xrange2[1]) then begin
;             djs_oplot, wave, smooth(flux,4), color=colors[ccount]
;             ccount++
;          endif
;       endfor
;       im_legend, [strtrim(spec.galaxy,2),'exp'+string(spec.frame,format='(I4.4)'),$
;         'z='+strtrim(string(spec.z,format='(F12.3)'),2)], $
;         /left, /top, box=0, charsize=1.7
;
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
;    im_plotconfig, psfile=qafile, /psclose, /gzip

return
end
