pro hizea_get_kband, out, clobber=clobber
; get K-band magnitudes for Alek's Keck run

    sdsspath = hizea_path(/sdss)
    isedpath = hizea_path(/ised)
    sfhgrid_basedir = hizea_path(/monte)
    sfhgrid_paramfile = hizea_path(/mass)+'hizea_sfhgrid.par'

    allsample = mrdfits(sdsspath+'hizea_galex_sdss_spitzer.fits.gz',1)
    keck = rsex(sdsspath+'keckao_sample.sex')

    ra = 15D*im_hms2dec(keck.ra) 
    dec = im_hms2dec(keck.dec)
    spherematch, allsample.ra, allsample.dec, ra, dec, 5.0/3600.0, m1, m2
    srt = sort(m2) & m1 = m1[srt] & m2 = m2[srt]
    sample = allsample[m1]
    ngal = n_elements(sample)

    out = struct_trimtags(sample,select=['ra','dec','z','galaxy','qso'])
    out = struct_addtags(out,replicate({kmag_kcorr: 0.0, kmag_kcorr_noirac: 0.0, $
      kmag_ised: 0.0, kmag_ised_noirac: 0.0, chi2_kcorr: 0.0, chi2_kcorr_noirac: 0.0, $
      chi2_ised: 0.0, chi2_ised_noirac: 0.0},ngal))
    
    djs_plot, ra, dec, psym=6, sym=2, xsty=3, ysty=3
    djs_oplot, allsample.ra, allsample.dec, psym=7, color='red'

; k-correct    
    kcorr = mrdfits(isedpath+'hizea_kcorrect.fits.gz',1,rows=m1)
    kcorr_noirac = mrdfits(isedpath+'hizea_kcorrect_noirac.fits.gz',1,rows=m1)

    out.kmag_kcorr = kcorr.k_mobs_fnuvugrizjhk[9]
    out.chi2_kcorr = kcorr.k_chi2
    out.kmag_kcorr_noirac = kcorr_noirac.k_mobs_fnuvugrizjhk[9]
    out.chi2_kcorr_noirac = kcorr_noirac.k_chi2
    
; isedfit    
    paramfile = isedpath+'hizea_isedfit.par'
    model = isedfit_restore(paramfile,ised,params=params,$
      iopath=isedpath,index=m1,sfhgrid_basedir=sfhgrid_basedir,/flam)
    model_noirac = isedfit_restore(paramfile,ised_noirac,params=params,$
      iopath=isedpath,index=m1,sfhgrid_basedir=sfhgrid_basedir,/flam,$
      outprefix='hizea_noirac')

    out.chi2_ised = ised.chi2
    out.chi2_ised_noirac = ised_noirac.chi2
    
    for ii = 0, ngal-1 do begin
       kmag = k_project_filters(k_lambda_to_edges(model[ii].wave),$
         model[ii].flux,filterlist='twomass_Ks.par')
       kmag_noirac = k_project_filters(k_lambda_to_edges(model_noirac[ii].wave),$
         model_noirac[ii].flux,filterlist='twomass_Ks.par')
       out[ii].kmag_ised = -2.5*alog10(reform(kmag))
       out[ii].kmag_ised_noirac = -2.5*alog10(reform(kmag_noirac))
    endfor

; write out
    outfile = sdsspath+'hizea_keckao_kband.fits'
    im_mwrfits, out, outfile, clobber=clobber
    
return
end
    
