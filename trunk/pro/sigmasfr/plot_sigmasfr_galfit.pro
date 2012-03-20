pro plot_sigmasfr_galfit, file, xcen, ycen, posa=posa, posb=posb, posc=posc, npix=npix
;+
;  amd120320 -- written to plot GALFIT data, model, and residuals
;
;-

if n_elements(xcen) lt 1L then xcen = 79.
if n_elements(ycen) lt 1L then ycen = 73.

if n_elements(posa) lt 4L then posa = [0.15, 0.6, 0.3, 0.8]
if n_elements(posb) lt 4L then posb = [0.3, 0.6, 0.45, 0.8]
if n_elements(posc) lt 4L then posc = [0.45, 0.6, 0.6, 0.8]

if n_elements(npix) lt 1L then npix = 50


data = readfits(file, header, ext=1)
model = readfits(file, header, ext=2)
resid = readfits(file, header, ext=3)

  stamp_data = dindgen(npix*2+1L,npix*2+1L)
  stamp_model = dindgen(npix*2+1L,npix*2+1L)
  stamp_resid = dindgen(npix*2+1L,npix*2+1L)

  for j = -1*npix, npix do begin
    for k = -1*npix, npix do begin
      stamp_data[j + npix,k + npix] = data[xcen + j,ycen + k]
      if stamp_data[j + npix,k + npix] ne $
        stamp_data[j + npix,k + npix] * 1. then $
        stamp_data[j + npix,k + npix] = 0.
      stamp_model[j + npix,k + npix] = model[xcen + j,ycen + k]
      if stamp_model[j + npix,k + npix] ne $
        stamp_model[j + npix,k + npix] * 1. then $
        stamp_model[j + npix,k + npix] = 0.
      stamp_resid[j + npix,k + npix] = resid[xcen + j,ycen + k]
      if stamp_resid[j + npix,k + npix] ne $
        stamp_resid[j + npix,k + npix] * 1. then $
        stamp_resid[j + npix,k + npix] = 0.
    endfor
  endfor


  ;dfpsplot, 'hst.ps', /landscape

  djs_iterstat, stamp_data, median=md, sigma=sig, sigrej=5.0
  min = md-1.*sig
  max = md+6.*sig
 
  imdisp, imgscl(stamp_data,max = max, min=min), $
    /usepos, /negative, bottom=40., pos=posa
  imdisp, imgscl(stamp_model,max = max, min=min), $
    /usepos, /negative, bottom=40., pos=posb
  imdisp, imgscl(stamp_resid,max = max, min=min), $
    /usepos, /negative, bottom=40., pos=posc

  ;dfpsclose

;stop

return
end
  
