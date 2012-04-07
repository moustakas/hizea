pro sigmasfr_lir, clobber=clobber
; jm12mar19ucsd - compute L(IR) for the SIGMASFR sample

    path = sigmasfr_path()

    cat = mrdfits(path+'sigmasfr_photometry.fits.gz',1)
    ngal = n_elements(cat)
    out = struct_addtags(struct_trimtags(cat,select=['galaxy','z']),$
      replicate({$
      lir_chary: -999.0, lir_dale: -999.0, lir_rieke: -999.0, $
      sfr_dale: -999.0, sfr_chary: -999.0, sfr_rieke: -999.0, $
      indx_chary: -999.0, indx_dale: -999.0, indx_rieke: -999.0},ngal))
    
    filters = strtrim(sigmasfr_filterlist(),2)
    these = where(filters eq 'wise_w3.par' or filters eq 'wise_w4.par')

    good = where(total(cat.maggies[these] gt 0,1) gt 0)
    out[good].lir_chary = alog10(im_wise2lir(cat[good].z,cat[good].maggies[these],$
      cat[good].ivarmaggies[these],/chary,model_indx=indx_chary))
    out[good].lir_dale = alog10(im_wise2lir(cat[good].z,cat[good].maggies[these],$
      cat[good].ivarmaggies[these],/dale,model_indx=indx_dale))
    out[good].lir_rieke = alog10(im_wise2lir(cat[good].z,cat[good].maggies[these],$
      cat[good].ivarmaggies[these],/rieke,model_indx=indx_rieke))

    out[good].indx_chary = indx_chary
    out[good].indx_dale = indx_dale
    out[good].indx_rieke = indx_rieke

    fact = 4.5D-44*im_lsun()
    out[good].sfr_chary = 10.0^out[good].lir_chary*fact
    out[good].sfr_dale = 10.0^out[good].lir_dale*fact
    out[good].sfr_rieke = 10.0^out[good].lir_rieke*fact

    im_mwrfits, out, path+'sigmasfr_lir.fits', clobber=clobber

    for i=0L, n_elements(good)-1L do begin
      SFR = wr11_lir_to_sfr(wr11_f24_to_lir($
        flux=cat[good[i]].maggies[10]*3631.d3, $
        zobj=cat[good[i]].z))
      splog, sfr, out[good[i]].sfr_chary
    endfor

return
end
    
