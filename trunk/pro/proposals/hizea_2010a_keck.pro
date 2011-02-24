function cmodel, x, p ; for METHOD 3
return, x # p
end

pro hizea_2010a_keck
; jm09sep04ucsd - make a pretty plot for the Keck/2010A proposal

; read the data    
    hizeapath = getenv('HIZEA_HIRES_DIR')+'/'
    fits = hizeapath+[$
      '09aug12/spec1d/J1558+3957_0238',$
      '09aug13/spec1d/J2140+1209_1195']+'.fits.gz'

; make the plot    
    psfile = getenv('PAPERSPATH')+'/proposals/Keck/2010A/hizea_hires_09aug.ps'
    im_plotconfig, 9, pos, psfile=psfile, xmargin=[1.2,0.3], charsize=1.9

    xrange1 = [2770,2812]
    yrange1 = [-0.09,1.4999]
    pad = 10
    colors = ['orange','blue','red','green','purple','cyan','yellow']
    this_order = [91,73]

    npoly = 6
    parinfo = {value: 1.0D}
    parinfo = replicate(parinfo,npoly)
    
    for ii = 0, 1 do begin
       splog, 'Reading '+fits[ii]
       spec = mrdfits(fits[ii],1,/silent)
       zobj = spec.z
; pick out the order we care about
       ww = (where(spec.order eq this_order[ii]))[0]
       good = where(spec.wave[*,ww] gt 0.0,ngood)
       good = good[pad:ngood-1-pad]
       npix = n_elements(good)
       wave = spec.wave[good,ww]/(1.0+zobj)
       flux = spec.flux[good,ww]
       ferr = spec.ferr[good,ww]
       case ii of
          0: mask = ((wave ge 2783) and (wave le 2803)) ne 1
          1: mask = ((wave ge 2787) and (wave le 2807)) ne 1
       endcase
       ivar = 1.0/ferr^2 * mask
       mgii = where(mask eq 0)
       splog, 'Median S/N = ', median(flux[mgii]/ferr[mgii])
;      print, minmax(wave)
; normalize
       pflux = poly_array(npix,npoly)
       coeff = mpfitfun('cmodel',pflux,flux,weights=ivar,$
         /quiet,parinfo=parinfo)
       cflux = pflux#coeff
       nflux = flux/cflux
;      djs_plot, wave, flux, yrange=[0,0.08]
;      djs_oplot, wave, cflux, color='blue'
;      cc = get_kbrd(1)
; make the plot
       if (ii eq 0) then begin
          xtickname1 = replicate(' ',10)
          xtitle = ''
       endif else begin
          delvarx, xtickname1
          xtitle1 = 'Rest Wavelength (\AA)'
       endelse
       
       djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=xrange1, $
         yrange=yrange1, xtickname=xtickname1, xtitle=xtitle1, $
         ytitle='', position=pos[*,ii], noerase=(ii ne 0), ytickinterval=0.5
       djs_oplot, wave, nflux, psym=10, color='dark grey'
;      djs_oplot, wave, smooth(nflux,2), psym=10
       djs_oplot, 2796.0*[1,1], !y.crange, line=5, color='blue', thick=8.0
       djs_oplot, 2803.0*[1,1], !y.crange, line=5, color='blue', thick=8.0
;      djs_oplot, !x.crange, [0,0], line=0, color='grey', thick=4.0

       im_legend, [strtrim(spec.galaxy,2),$
         'z = '+strtrim(string(spec.z,format='(F12.3)'),2)], $
;        't_{exp} = '+string(spec.exp,format='(I4)')+' s'], $
         /left, /bottom, box=0, charsize=1.5
    endfor

    xyouts, pos[0,0]-0.07, pos[1,0], 'Relative Flux', align=0.5, $
      orientation=90, charsize=!p.charsize, /norm
    
    im_plotconfig, psfile=qafile, /psclose, /gzip

return
end
    
