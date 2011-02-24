function irac_phot, image, nocntrd=nocntrd
;+
;
; amd101208 -- modified from previous irac photometry code.  here the
;              goal is to measure 3.6 and 4.5 micron photometry for
;              100 post-starburst galaxies.
;
; amd110224 -- added this comment to see if i can commit via the google 
;              svn repository
;-            



im = readfits(image, header, /silent)

object = sxpar(header, 'object')
ra_hms = strmid(object, 5,2)+':'+strmid(object, 7,2)+':'+strmid(object, 9,5)
dec_hms = strmid(object, 14,3)+':'+strmid(object, 17,2)+':'+strmid(object, 19,4)

extast, header, astr
ad2xy, im_hms2dec(ra_hms)*15., im_hms2dec(dec_hms), astr, xcen, ycen
  
cntrd, im, xcen, ycen, x, y, 4.
xin = xcen
yin = ycen
xout = x
yout = y
xy2ad, xin, yin, astr, rain, decin
xy2ad, xout, yout, astr, raout, decout

if not keyword_set(nocntrd) then begin
  xcen = x
  ycen = y
endif


fluxconv = sxpar(header, 'fluxconv') 
;exptime = sxpar(header, 'exptime')
exptime=1
gain = sxpar(header, 'gain') 
factor = (!pi/180.d)^2.d * (1.d/3600.d)^2.d * $
  (0.6d / 1.d)^2.d * 1.d6 * 1.d-23 ; MJy/sr -> erg/s/cm^-2/Hz/pix

im = im/fluxconv*exptime*gain

; to get from units of image data (MJy/sr) to electrons, 
; (1) divide by FLUXCONV (BUNIT to DN/s)
; (2) multiply by EXPTIME (DN/s to DN)
; (3) mulitply by GAIN (DN to e)

rd = findgen(20.)+1
badpix = [-32767., 80000.d]
skyrad = [24.,40.]
aper, im, xcen, ycen, flux, error, sky, skyerr, gain, rd, $
  skyrad, badpix, /flux, /silent

if flux[0] ne flux[0]*1. then begin
  splog, object, ' needs to chill out ...'
  ;splog, ra_hms, dec_hms
  return, -1L
endif

flux = flux/gain/exptime*fluxconv*factor*1.d29
error = error/gain/exptime*fluxconv*factor*1.d29 ; convert to uJy
;flux = flux/exptime*fluxconv*factor*1.d26
;error = error/exptime*fluxconv*factor*1.d26


posa = [0.15,0.1,0.9,0.45]
posb = [0.525,0.6,0.9,0.9]
posc = [0.1,0.7,0.475,0.9]
color = ['white', 'orange', 'red', 'blue', 'dark green', 'brown', 'cyan', 'black', 'black']

plot, rd, flux, $
  xtitle='Aperture radius [pixels]', $
  xrange=[0.,max(rd)], xsty=1, $
  ytitle=textoidl('f_{\nu} [\muJy]'), $
  pos=posa, charsize=1.5, psym=3
oploterror, rd, flux, error, psym=3



  xyouts, 0.10, 0.82, $
    '(rain , decin ) =('+string(im_dec2hms(rain/15., /colon), format='(a11)')+$
    ', '+string(im_dec2hms(decin, /colon), format='(a11)')+')', /normal
  xyouts, 0.10, 0.80, $
    '(raout, decout)=('+string(im_dec2hms(raout/15., /colon), format='(a11)')+$
    ', '+string(im_dec2hms(decout, /colon), format='(a11)')+')', /normal

  deez = amd_radec_dist(ra1_deg=rain, dec1_deg=decin, $
    ra2_deg=raout, dec2_deg=decout)
  xyouts, 0.10, 0.78, $
    'dif (arcsec): '+string(deez.dist_arc, format='(f5.2)'), /normal

if keyword_set(nocntrd) then str='<-' else str=''

  xyouts, 0.10, 0.74, $
    '(xin , yin )=('+string(xin, format='(f8.3)')+', '+$
    string(yin, format='(f8.3)')+')'+str, /normal
 
if keyword_set(nocntrd) then str='' else str='<-'

  xyouts, 0.10, 0.72, $
    '(xout, yout)=('+string(xout, format='(f8.3)')+', '+$
    string(yout, format='(f8.3)')+')'+str, /normal
  xyouts, 0.10, 0.70, $
    'dif (x, y): '+'('+string(xout-xin, format='(f8.3)')+$
    string(yout-yin, format='(f8.3)')+')', /normal





rel = where(rd eq 4., nrel)
fpoint = flux[rel]*1.213
plots, !x.crange, [fpoint, fpoint]
xyouts, 0.15, 0.46, $
  'fpoint='+string(fpoint, format='(f8.1)')+$
  textoidl('\pm')+string(error[rel]*1.213, format='(f6.1)'), /normal


rel = where(rd eq 6., nrel)
fpoint = flux[rel]*1.1125
plots, !x.crange, [fpoint, fpoint]
xyouts, 0.15, 0.48, $
  'fpoint='+string(fpoint, format='(f8.1)')+$
  textoidl('\pm')+string(error[rel]*1.1125, format='(f6.1)'), /normal

rel = where(rd eq 10., nrel)
fpoint = flux[rel]
plots, !x.crange, [fpoint, fpoint]
xyouts, 0.15, 0.50, $
  'fpoint='+string(fpoint, format='(f8.1)')+$
  textoidl('\pm')+string(error[rel], format='(f6.1)'), /normal



xyouts, 0.1,0.92, object, /normal, size=2.0
directory = strmid(image,0,9) 
splog, directory
xyouts, 0.1, 0.88, directory, /normal, size=1.5

;0.15, 0.46

zody_est = sxpar(header, 'zody_est') 
xyouts, 0.1, 0.60, 'sky: '+string(sky/exptime*fluxconv, format='(f6.3)')+$
  textoidl('\pm')+$
  string(skyerr/exptime*fluxconv, format='(f6.3)')+$
  ' ('+string(zody_est, format='(f6.3)')+')', /normal


  npix = 40.
  stamp = dindgen(npix*2+1L,npix*2+1L)
  for j = -1*npix, npix do begin
    for k = -1*npix, npix do begin
      stamp[k + npix,j + npix] = im[xcen + k,ycen + j]
      if stamp[k + npix,j + npix] ne stamp[k + npix,j + npix] * 1. then $
        stamp[k + npix,j + npix] = 0.
    endfor
  endfor



  djs_iterstat, stamp, median=md, sigma=sig, sigrej=5.0
  ; imgscl, logscl, hist_equal
  imdisp, imgscl(stamp,max = md+15.*sig, min=md+0.*sig), $
    pos=posb, /usepos, /negative, bottom=40.

  plot, [0.,0.], [0.,0.], pos=posb, /noerase, $
    xrange=[xcen-npix,xcen+npix], xsty=5, $
    yrange=[ycen-npix,ycen+npix], ysty=5

meh = [4.,6.,10.,20.,24.,40.]

  for i=0L, n_elements(meh)-1L do begin
    tvcircle, meh[i], xcen, ycen, /data, color=djs_icolor(color[i])
  endfor



  tvbox, npix*2+1L, xcen, ycen, /data, color=djs_icolor('black')

;dfpsclose

;stop 

rel = where(rd eq 4., nrel)
fout = flux[rel]*1.213


return, fout
end
