;+
;  amd120207 -- written to creates maggies and ivarmaggies input for
;               input to ISEDFIT
;
;  
;
;-

nobj = 1.
nbands = 9.

range = [4000., 8000.]
nbins = 8.
dx = (range[1] - range[0])/nbins


output = {                     $
  ra:   0.d,                   $
  dec:  0.d,                   $
  z:    0.d,                   $
  galaxy: ' ',                 $
  filterlist: strarr(nbands+nbins),  $
  maggies: dblarr(nbands+nbins),     $ 
  ivarmaggies: dblarr(nbands+nbins)  }

output.ra = 136.3483394d
output.dec = 57.98680093d
output.z = 0.7116d

output.filterlist[0:nbands-1L] = ['galex_FUV.par','galex_NUV.par','sdss_u0.par','sdss_g0.par','sdss_r0.par','sdss_i0.par','sdss_z0.par','spitzer_irac_ch1.par','spitzer_irac_ch2.par']

output.galaxy = hogg_iau_name(output.ra,output.dec,'SDSS')

phot = rsex('background_photometry.dat')

mags = -2.5d*alog10(phot.fnu/3631.d-23)
maggies = 10.^(mags*(-1.)/2.5)
output.maggies[0:nbands-1L] = maggies

sigmag = -2.5d*alog10(phot.error/3631.d-23)
sigmaggies = 10.^(sigmag*(-1.)/2.5)
ivarmaggies = 1. / sigmaggies^2.
output.ivarmaggies[0:nbands-1L] = ivarmaggies

yo = mrdfits('j0905+5759_lris_flux_v120125.fits', 1)

dfpsplot, 'spectrophotometry.ps', /landscape, /color

  plot, yo.wavelength, yo.flux, psym=10
  c = 3.d8 * 1.d10
  fnu = yo.flux*1.d-17 * (yo.wavelength)^2 / c
  sigfnu = yo.error*1.d-17 * (yo.wavelength)^2 / c
  plot, yo.wavelength, fnu, psym=10
  oplot, phot.wave, phot.fnu, psym=6, symsize=1, $
    color=djs_icolor('red')



  wave = fltarr(nbins)
  flux = fltarr(nbins)

  for i=0L, nbins-1L do begin
    rel = where(yo.wavelength gt range[0]+dx*i-250. and $
                yo.wavelength lt range[0]+dx*(i+1L)+250., nrel)
    wave[i] = range[0]+dx*i + 225.
    ;meanerr, fnu[rel], sigfnu[rel], mean, sigmam, sigmad 
    ;flux[i] = mean
    p = [range[0]+dx*i + 225., 250., 10.]
    yvals = gauss1(yo[rel].wavelength, p)
    splog, total(yvals)
    oplot, yo[rel].wavelength, yvals * median(phot.fnu) / median(yvals)
    openw, 1, 'amd_'+string(wave[i], format='(i4)')+'_v2.par'
    printf, 1, 'typedef struct {'
    printf, 1, ' float lambda;'
    printf, 1, ' float pass;'
    printf, 1, '} KFILTER;'
    printf, 1, '' 
    ;flux[i] = 0.
    for j=0L, nrel-1L do begin
      printf, 1, 'KFILTER '+string(yo[rel[j]].wavelength, format='(f7.2)')+$
        ' '+string(yvals[j], format='(f6.4)')
      flux[i] = flux[i] + fnu[rel[j]] * yvals[j] / total(yvals)
     endfor
    close, 1
  endfor

;; 4000-4500
;  b4250 = where(yo.wavelength gt 4000. and yo.wavelength lt 4500., nb4250)
;  meanerr, fnu[b4250], sigfnu[b4250], f4250, sigmam, sigmad
;; 4500-5000
;  b4750 = where(yo.wavelength gt 4500. and yo.wavelength lt 5000., nb4750)
;  meanerr, fnu[b4750], sigfnu[b4750], f4750, sigmam, sigmad
;
;  wave = [4250., 4750.]
;  flux = [f4250, f4750]

  oplot, wave, flux,  psym=6, color=djs_icolor('blue') 

  spawn, 'tar cvzf amd_filters_v2.tar.gz amd*v2.par'

;filters =  file_search('amd*.par')
output.filterlist[nbands:nbands+nbins-1L] = file_search('amd*v2.par')

mags = -2.5d*alog10(flux/3631.d-23)
maggies = 10.^(mags*(-1.)/2.5)
output.maggies[nbands:nbands+nbins-1L] = maggies

sigmag = -2.5d*alog10(flux/60./3631.d-23)
sigmaggies = 10.^(sigmag*(-1.)/2.5)
ivarmaggies = 1. / sigmaggies^2.
output.ivarmaggies[nbands:nbands+nbins-1L] = ivarmaggies

mwrfits, output, 'j0905+5759_isedfit_input_v2.fits.gz'

dfpsclose

stop

end
