pro build_hizea_photometry, out, clobber=clobber
; jm11apr06ucsd - build the combined GALEX+SDSS+Spitzer photometric
; dataset for the full hizea sample

    hizeapath = getenv('HIZEA_DATA')+'/'
    outfile = hizeapath+'sdss/hizea_galex_sdss_spitzer.fits'

    info = mrdfits(hizeapath+'sdss/hizea_simplefit_dr7.fit',1)
    sdss = mrdfits(hizeapath+'sdss/hizea_photo_dr7.fit',1)
    galex = mrdfits(hizeapath+'sdss/hizea_galex_gr6.fits.gz',1)
    ngal = n_elements(info)

; line-match everything
    irac1 = mrdfits(hizeapath+'spitzer/hizea_spitzer.fits.gz',1)
    irac = im_empty_structure(irac1,empty_value=-999.0,ncopies=ngal)
    spherematch, info.ra, info.dec, irac1.ra, irac1.dec, 1.0/3600.0, m1, m2
    srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
    irac[m1] = irac1[m2]

; pack it all together
    out = struct_addtags(struct_trimtags(info,select=['ra','dec','z']),$
      replicate({galaxy: '', maggies: fltarr(9), ivarmaggies: fltarr(9), $
      qso: 0},ngal))
    out.galaxy = hogg_iau_name(out.ra,out.dec,'SDSS')

; is this in the DR7 quasar catalog?    
    qq = mrdfits(getenv('CATALOGS_DIR')+'/10schneider/10schneider.fits.gz',1)
    spherematch, qq.raj2000, qq.dej2000, out.ra, out.dec, 1.0/3600.0, m1, m2
    out[m2].qso = 1
    
; galex
    im_galex_to_maggies, galex, maggies, ivarmaggies
    out.maggies[0:1] = maggies
    out.ivarmaggies[0:1] = ivarmaggies

; sdss    
    sdss_to_maggies, maggies, ivarmaggies, tsobj=sdss, flux='model'
    out.maggies[2:6] = maggies
    out.ivarmaggies[2:6] = ivarmaggies

; irac; the IRAC fluxes are in microJy, so to convert to maggies
; multiply by 10^(-0.4*23.9)
    imaggies = fltarr(2,ngal)
    iivarmaggies = fltarr(2,ngal)

    ch1good = where(irac.ch1flux gt 0.0)
    ch2good = where(irac.ch2flux gt 0.0)
    imaggies[0,ch1good] = irac[ch1good].ch1flux*10.0^(-0.4*23.9)
    iivarmaggies[0,ch1good] = 1.0/(imaggies[0,ch1good]*0.1)^2
    
    imaggies[1,ch2good] = irac[ch2good].ch2flux*10.0^(-0.4*23.9)
    iivarmaggies[1,ch2good] = 1.0/(imaggies[1,ch2good]*0.1)^2

    out.maggies[7:8] = imaggies
    out.ivarmaggies[7:8] = iivarmaggies

    im_mwrfits, out, outfile, clobber=clobber
    
return
end
