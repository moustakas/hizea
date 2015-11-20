pro build_hizea_photometry, out, clobber=clobber
; jm11apr06ucsd - build the combined GALEX+SDSS+Spitzer photometric
;   dataset for the full hizea sample
; jm15nov16siena - updated to the latest sample and photometric data model

    common com_phot, sdss, unwise

    photpath = getenv('HIZEA_DATA')+'/phot/'
    hizeadir = massprofiles_path(/code)+'etc/'
    
    filt = hizea_filterlist()
    weff = k_lambda_eff(filterlist=filt)
    nfilt = n_elements(filt)

    sample = rsex(hizeadir+'hizea_sample2.txt')
    ngal = n_elements(sample)

    phot = struct_trimtags(sample,select=['GALAXY','RA','DEC','Z'])
    phot = struct_addtags(phot,replicate({maggies: fltarr(nfilt), $
      ivarmaggies: fltarr(nfilt)},ngal))
    
; galex
    galex = mrdfits(photpath+'hizea_sample2_galex_gr6.fits',1)
    im_galex_to_maggies, galex, maggies, ivarmaggies
    phot.maggies[0:1] = maggies
    phot.ivarmaggies[0:1] = ivarmaggies
    
;; spitzer
;    spit = mrdfits(photpath+'hizea_spitzer.fits.gz',1)
;    spherematch, sample.ra, sample.dec, spit.ra, spit.dec, 1D/3600, m1, m2
;
;    spit = rsex(hizeadir+'hst_sample_spitzer.dat')
;    gal = strmid(phot.galaxy,0,5)+strmid(phot.galaxy,10,1)+$
;      strmid(phot.galaxy,124)
;    match, strtrim(spit.galaxy,2), strmid(gal,2), m1, m2
;    ww = where(spit[m1].ch1_flux gt 0)
;    phot[m2[ww]].maggies[7:8] = transpose([[spit[m1[ww]].ch1_flux],$
;      [spit[m1[ww]].ch2_flux]]*10^(-0.4*23.9))
;    phot[m2[ww]].ivarmaggies[7:8] = 1D/(transpose([[spit[m1[ww]].ch1_sig],$
;      [spit[m1[ww]].ch2_sig]]*10^(-0.4*23.9)))^2
    
;; galex; take the weighted mean of the GALEX NUV observations
;    galex = im_read_tbl(path+'sigmasfr_schiminovich_galex.txt',/silent)
;    im_galex_to_maggies, galex, mm, ii, /psf
;
;    nuv = im_weighted_mean(mm[1,*],error=1D/sqrt(ii[1,*]),wsigma=nuverr)
;    phot.maggies[0] = reform(mm[0,3])
;    phot.ivarmaggies[0] = reform(ii[0,3])
;    phot.maggies[1] = nuv
;    phot.ivarmaggies[1] = 1D/nuverr^2

;; galex (pre-schiminovich)
;    galex = im_read_tbl(sigmasfrdir+'galexview_output.tbl',/silent)
;    im_galex_to_maggies, galex, mm, ii, /psf
;; loop through the GALEX table
;    for jj=0L, n_elements(phot)-1L do begin
;       match=where(galex.UPLOADSEARCHID eq jj+1L, nmatch)
;; find useful FUV data
;       fuvgood=where(mm[0,match] gt 0., nfuvgood)
;       if nfuvgood lt 1L then begin
;          phot[jj].maggies[0] = 0.
;          phot[jj].ivarmaggies[0] = 0.
;       endif else begin
;          if nfuvgood eq 1L then begin
;             phot[jj].maggies[0] = mm[0,match[fuvgood]]
;             phot[jj].ivarmaggies[0] = ii[0,match[fuvgood]]      
;          endif else begin
;             fuv = im_weighted_mean(mm[0,match[fuvgood]],$
;               error=1D/sqrt(ii[0,match[fuvgood]]),wsigma=fuverr)
;             phot[jj].maggies[0] = fuv
;             phot[jj].ivarmaggies[0] = 1D/fuverr^2
;          endelse
;       endelse
;; find useful NUV data
;       nuvgood=where(mm[1,match] gt 0., nnuvgood)
;       if nnuvgood lt 1L then begin
;          phot[jj].maggies[1] = 0.
;          phot[jj].ivarmaggies[1] = 0.
;       endif else begin
;          if nnuvgood eq 1L then begin
;             phot[jj].maggies[1] = mm[1,match[nuvgood]]
;             phot[jj].ivarmaggies[1] = ii[1,match[nuvgood]]      
;          endif else begin
;             nuv = im_weighted_mean(mm[1,match[nuvgood]],$
;               error=1D/sqrt(ii[1,match[nuvgood]]),wsigma=nuverr)
;             phot[jj].maggies[1] = nuv
;             phot[jj].ivarmaggies[1] = 1D/nuverr^2
;          endelse
;       endelse  
;    endfor

;; is this in the DR7 quasar catalog?    
;    qq = mrdfits(getenv('CATALOGS_DIR')+'/10schneider/10schneider.fits.gz',1)
;    spherematch, qq.raj2000, qq.dej2000, out.ra, out.dec, 1.0/3600.0, m1, m2
;    out[m2].qso = 1

; sdss/dr12
    if (n_elements(sdss) eq 0L) then sdss = $
      mrdfits(getenv('IM_ARCHIVE_DIR')+'/data/sdss/dr12/photoPosPlate-dr12.fits',1,$
      columns=['ra','dec'])
    spherematch, sdss.ra, sdss.dec, phot.ra, phot.dec, 1D/3600, m1, m2
    sdssphot = mrdfits(getenv('IM_ARCHIVE_DIR')+'/data/sdss/dr12/photoPosPlate-dr12.fits',1,rows=m1)
    sdss_to_maggies, mm, ii, calib=sdssphot, flux='model'

    phot[m2].maggies[2:6] = mm
    phot[m2].ivarmaggies[2:6] = ii

; unwise
    if (n_elements(unwise) eq 0L) then unwise = $
      mrdfits(getenv('IM_ARCHIVE_DIR')+'/data/sdss/dr10/specmatch-dr10.fits',1)
    spherematch, unwise.plug_ra, unwise.plug_dec, phot.ra, phot.dec, 5D/3600, m1, m2
    unwise_to_maggies, unwise[m1], mm, ii, prefix='wise', ratag='plug_ra', $
      dectag='plug_dec'
    
    phot[m2].maggies[9:12] = mm
    phot[m2].ivarmaggies[9:12] = ii

; write out
;   minerr = replicate(0.05,nfilt)
;   maggies = phot.maggies
;   ivarmaggies = phot.ivarmaggies
;   k_minerror, maggies, ivarmaggies, minerr
;   phot.maggies = maggies
;   phot.ivarmaggies = ivarmaggies
    im_mwrfits, phot, photpath+'hizea_sample2_phot.fits', /clob

stop
    
return
end
    
    
;    hizeapath = getenv('HIZEA_DATA')+'/'
;    outfile = hizeapath+'sdss/hizea_galex_sdss_spitzer.fits'
;
;    info = mrdfits(hizeapath+'sdss/hizea_simplefit_dr7.fit',1)
;    sdss = mrdfits(hizeapath+'sdss/hizea_photo_dr7.fit',1)
;    galex = mrdfits(hizeapath+'sdss/hizea_galex_gr6.fits.gz',1)
;    ngal = n_elements(info)
;
;; line-match everything
;    irac1 = mrdfits(hizeapath+'spitzer/hizea_spitzer.fits.gz',1)
;    irac = im_empty_structure(irac1,empty_value=-999.0,ncopies=ngal)
;    spherematch, info.ra, info.dec, irac1.ra, irac1.dec, 1.0/3600.0, m1, m2
;    srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
;    irac[m1] = irac1[m2]
;
;; pack it all together
;    out = struct_addtags(struct_trimtags(info,select=['ra','dec','z']),$
;      replicate({galaxy: '', maggies: fltarr(9), ivarmaggies: fltarr(9), $
;      qso: 0},ngal))
;    out.galaxy = hogg_iau_name(out.ra,out.dec,'SDSS')
;
;    
;; galex
;    im_galex_to_maggies, galex, maggies, ivarmaggies
;    out.maggies[0:1] = maggies
;    out.ivarmaggies[0:1] = ivarmaggies
;
;; sdss    
;    sdss_to_maggies, maggies, ivarmaggies, tsobj=sdss, flux='model'
;    out.maggies[2:6] = maggies
;    out.ivarmaggies[2:6] = ivarmaggies
;
;; irac; the IRAC fluxes are in microJy, so to convert to maggies
;; multiply by 10^(-0.4*23.9)
;    imaggies = fltarr(2,ngal)
;    iivarmaggies = fltarr(2,ngal)
;
;    ch1good = where(irac.ch1flux gt 0.0)
;    ch2good = where(irac.ch2flux gt 0.0)
;    imaggies[0,ch1good] = irac[ch1good].ch1flux*10.0^(-0.4*23.9)
;    iivarmaggies[0,ch1good] = 1.0/(imaggies[0,ch1good]*0.1)^2
;    
;    imaggies[1,ch2good] = irac[ch2good].ch2flux*10.0^(-0.4*23.9)
;    iivarmaggies[1,ch2good] = 1.0/(imaggies[1,ch2good]*0.1)^2
;
;    out.maggies[7:8] = imaggies
;    out.ivarmaggies[7:8] = iivarmaggies
;
;    im_mwrfits, out, outfile, clobber=clobber, /nogzip
;    
;return
;end
