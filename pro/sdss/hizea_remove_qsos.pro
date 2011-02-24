pro hizea_remove_qsos, clobber=clobber
; jm10dec22ucsd - remove QSOs from the parent HIZEA sample

    sdsspath = getenv('HIZEA_DATA')+'/sdss/'
    phot = mrdfits(sdsspath+'hizea_photo_dr7.fit',1)
    fit = mrdfits(sdsspath+'hizea_simplefit_dr7.fit',1)
    
    qq = mrdfits(getenv('CATALOGS_DIR')+'/10schneider/10schneider.fits.gz',1)
    spherematch, qq.raj2000, qq.dej2000, phot.ra, phot.dec, 1.0/3600.0, m1, m2

    keep = lindgen(n_elements(phot))
    remove, m2, keep
    outphot = phot[keep]
    outfit = fit[keep]

    srt = sort(outphot.ra)
    outphot = outphot[srt]
    outfit = outfit[srt]

; special case
    keep = where(outphot.ra ne 0.0)
    outphot = outphot[keep]
    outfit = outfit[keep]
    
    im_mwrfits, outphot, sdsspath+'hizea_sdss_photo_dr7.fits', clobber=clobber
    im_mwrfits, outfit, sdsspath+'hizea_sdss_simplefit_dr7.fits', clobber=clobber
    
return
end
    
