pro plot_sigmasfr_galfit, file, xcen, ycen, posa=posa, posb=posb, posc=posc, npix=npix, onlydata=onlydata
;+
;  amd120320 -- written to plot GALFIT data, model, and residuals
;
;-

if n_elements(xcen) lt 1L then xcen = 79.
if n_elements(ycen) lt 1L then ycen = 73.

; fiducial spacing on the top left
;if n_elements(posa) lt 4L then posa = [0.15, 0.65, 0.35, 0.9]
;if n_elements(posb) lt 4L then posb = [0.35, 0.65, 0.55, 0.9]
;if n_elements(posc) lt 4L then posc = [0.55, 0.65, 0.75, 0.9]

; squeeze things horizontally on the top left
;if n_elements(posa) lt 4L then posa = [0.1, 0.7, 0.3, 0.95]
;if n_elements(posb) lt 4L then posb = [0.3, 0.7, 0.5, 0.95]
;if n_elements(posc) lt 4L then posc = [0.5, 0.7, 0.7, 0.95]

; squeeze and compress horizontally on the top left
;if n_elements(posa) lt 4L then posa = [0.1, 0.75, 0.25, 0.95]
;if n_elements(posb) lt 4L then posb = [0.25, 0.75, 0.4, 0.95]
;if n_elements(posc) lt 4L then posc = [0.4, 0.75, 0.55, 0.95]

; vertical spacing on bottom right
;if n_elements(posa) lt 4L then posa = [0.8, 0.5, 0.95, 0.7]
;if n_elements(posb) lt 4L then posb = [0.8, 0.3, 0.95, 0.5]
;if n_elements(posc) lt 4L then posc = [0.8, 0.1, 0.95, 0.3]

; horizontal spacing on the bottom right
if n_elements(posa) lt 4L then posa = [0.26, 0.13, 0.49, 0.39]
if n_elements(posb) lt 4L then posb = [0.49, 0.13, 0.72, 0.39]
if n_elements(posc) lt 4L then posc = [0.72, 0.13, 0.95, 0.39]

if n_elements(npix) lt 1L then npix = 50


data = readfits(file, header, ext=1)
model = readfits(file, header, ext=2)
resid = readfits(file, header, ext=3)

gain = 1
rd = findgen(npix-10L)+1L
badpix = [-32767., 80000.d]
skyrad = [npix-10L, npix]
aper, data, xcen, ycen, fluxdata, error, sky, skyerr, gain, rd, $
  skyrad, badpix, /flux, /silent
aper, model, xcen, ycen, fluxmodel, error, sky, skyerr, gain, rd, $
  skyrad, badpix, /flux, /silent
aper, resid, xcen, ycen, fluxresid, error, sky, skyerr, gain, rd, $
  skyrad, badpix, /flux, /silent
splog, max(fluxmodel) / max(fluxdata), max(fluxresid) / max(fluxdata)
splog, fluxmodel[9] / fluxdata[9], fluxresid[9] / fluxdata[9]

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
  max = md+3*sig
  imdisp, imgscl(stamp_data,max = max, min=min), $
    /usepos, /negative, bottom=40., pos=posa

if not keyword_set(onlydata) then begin

  ;djs_iterstat, stamp_model, median=md, sigma=sig, sigrej=5.0
  ;min = md-1.*sig
  ;max = md+4.*sig
  imdisp, imgscl(stamp_model,max = max, min=min), $
    /usepos, /negative, bottom=40., pos=posb

  djs_iterstat, stamp_resid, median=md, sigma=sig, sigrej=5.0
  min = md-1.*sig
  max = md+4.*sig
  imdisp, imgscl(stamp_resid,max = max, min=min), $
    /usepos, /negative, bottom=40., pos=posc

endif

  ;dfpsclose

;stop

return
end
  
