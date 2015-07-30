pro massprofiles_lir, debug=debug, clobber=clobber
; jm15jul06siena - compute L(IR) for the MASSPROFILES sample

    path = massprofiles_path()

    outfile = path+'massprofiles_lir_all.fits'
    outsedsfile = path+'massprofiles_lir_seds_all.fits'
    
;   kcorr = mrdfits(path+'massprofiles_kcorrect.fits.gz',1)
    cat = mrdfits(path+'massprofiles_photometry.fits.gz',1)
    ngal = n_elements(cat)

; nuLnu is the rest-frame luminosity 
    out = struct_addtags(struct_trimtags(cat,select=['galaxy','z']),$
      replicate({$
      lir_chary: -999.0, lir_dale: -999.0, lir_rieke: -999.0, $
      nuLnu_chary: dblarr(3)-999.0, nuLnu_dale: dblarr(3)-999.0, nuLnu_rieke: dblarr(3)-999.0, $
      ivar_nulnu_chary: dblarr(3)-999.0, ivar_nulnu_dale: dblarr(3)-999.0, ivar_nulnu_rieke: dblarr(3)-999.0, $
      sfr_chary: -999.0, sfr_dale: -999.0, sfr_rieke: -999.0, $
;     uvsfr_chary: -999.0, uvsfr_dale: -999.0, uvsfr_rieke: -999.0, $
      indx_chary: -999.0, indx_dale: -999.0, indx_rieke: -999.0, $
      irx_chary: -999.0, irx_dale: -999.0, irx_rieke: -999.0, $
      a1500_chary: -999.0, a1500_dale: -999.0, a1500_rieke: -999.0},ngal))

    outseds = replicate({$
      modelwave_chary: fltarr(1366), modelwave_dale: fltarr(1930), modelwave_rieke: fltarr(6270), $
      modelmab_chary: fltarr(1366), modelmab_dale: fltarr(1930), modelmab_rieke: fltarr(6270)},$
      ngal)
    
    filters = strtrim(massprofiles_filterlist(/allbands),2)
    these = where(filters eq 'wise_w3.par' or filters eq 'wise_w4.par')

; throw out upper limits to make sure the K-corrections in
; im_simple_kcorrect get computed from the band where the object was
; detected 
    good = where(total(cat.maggies[these] gt 0,1) gt 0)
    wmaggies = cat[good].maggies[these]
;   wivarmaggies = cat[good].ivarmaggies[these]
    wivarmaggies = cat[good].ivarmaggies[these]*(wmaggies gt 0)
    
    out[good].lir_chary = alog10(im_wise2lir(cat[good].z,wmaggies,wivarmaggies,$
      /chary,model_indx=indx_chary,nulnu=nulnu_chary,debug=debug,$
      ivar_nulnu=ivar_nulnu_chary,modelwave=modelwave_chary,modelmab=modelmab_chary))
    out[good].lir_dale = alog10(im_wise2lir(cat[good].z,wmaggies,wivarmaggies,$
      /dale,model_indx=indx_dale,nulnu=nulnu_dale,$
      ivar_nulnu=ivar_nulnu_dale,modelwave=modelwave_dale,modelmab=modelmab_dale))
    out[good].lir_rieke = alog10(im_wise2lir(cat[good].z,wmaggies,wivarmaggies,$
      /rieke,model_indx=indx_rieke,nulnu=nulnu_rieke,$
      ivar_nulnu=ivar_nulnu_rieke,modelwave=modelwave_rieke,modelmab=modelmab_rieke))

    outseds[good].modelwave_chary = modelwave_chary
    outseds[good].modelmab_chary = modelmab_chary
    outseds[good].modelwave_dale = modelwave_dale
    outseds[good].modelmab_dale = modelmab_dale
    outseds[good].modelwave_rieke = modelwave_rieke
    outseds[good].modelmab_rieke = modelmab_rieke

    out[good].indx_chary = indx_chary
    out[good].indx_dale = indx_dale
    out[good].indx_rieke = indx_rieke

    out[good].nulnu_chary = nulnu_chary
    out[good].nulnu_dale = nulnu_dale
    out[good].nulnu_rieke = nulnu_rieke

    out[good].ivar_nulnu_chary = ivar_nulnu_chary
    out[good].ivar_nulnu_dale = ivar_nulnu_dale
    out[good].ivar_nulnu_rieke = ivar_nulnu_rieke

    fact = 4.5D-44*im_lsun()
    out[good].sfr_chary = 10.0^out[good].lir_chary*fact
    out[good].sfr_dale = 10.0^out[good].lir_dale*fact
    out[good].sfr_rieke = 10.0^out[good].lir_rieke*fact

;    for i=0L, n_elements(good)-1L do begin
;      SFR = wr11_lir_to_sfr(wr11_f24_to_lir($
;        flux=cat[good[i]].maggies[10]*3631.d3, $
;        zobj=cat[good[i]].z))
;      splog, sfr, out[good[i]].sfr_chary
;    endfor

;; compute upper limits using just the 22-micron flux (ignore the
;; 12-micron flux, if any)
;    crap = where((cat.maggies[these[1]] eq 0) and (cat.ivarmaggies[these[1]] gt 0))
;;   crap = where((total(cat.maggies[these] gt 0,1) eq 0) and (cat.ivarmaggies[these[1]] gt 0))
;    wmaggies = 1.0/sqrt(cat[crap].ivarmaggies[these])
;    wivarmaggies = 1.0/(0.05*wmaggies)^2
;    wivarmaggies[0,*] = 0
;
;    lir_chary = alog10(im_wise2lir(cat[crap].z,wmaggies,wivarmaggies,$
;      /chary,model_indx=indx_chary,nulnu=nulnu_chary,$
;      ivar_nulnu=ivar_nulnu_chary,modelwave=modelwave_chary,modelmab=modelmab_chary))
;    lir_dale = alog10(im_wise2lir(cat[crap].z,wmaggies,wivarmaggies,$
;      /dale,model_indx=indx_dale,nulnu=nulnu_dale,$
;      ivar_nulnu=ivar_nulnu_dale,modelwave=modelwave_dale,modelmab=modelmab_dale))
;    lir_rieke = alog10(im_wise2lir(cat[crap].z,wmaggies,wivarmaggies,$
;      /rieke,model_indx=indx_rieke,nulnu=nulnu_rieke,$
;      ivar_nulnu=ivar_nulnu_rieke,modelwave=modelwave_rieke,modelmab=modelmab_rieke))
;
;    out[crap].lir_chary = -lir_chary
;    out[crap].lir_dale = -lir_dale
;    out[crap].lir_rieke = -lir_rieke
;
;    out[crap].nulnu_chary = -nulnu_chary
;    out[crap].nulnu_dale = -nulnu_dale
;    out[crap].nulnu_rieke = -nulnu_rieke

;   niceprint, cat.galaxy, out.lir_chary, out.nulnu_chary, cat.maggies[these[1]], cat.ivarmaggies[these[1]]
    
; write out    
    im_mwrfits, out, outfile, clobber=clobber
    im_mwrfits, outseds, outsedsfile, clobber=clobber

; subselect the main sample
    these = 'J'+['1506+5402','0905+5759','1341-0321','0944+0930',$
      '2140+1209','0826+4305','1613+2834','1219+0336','1107+0417',$
      '0901+0314','0106-1023','2116-0634']
    match, out.galaxy, these, m1, m2
    srt = sort(m1)
    im_mwrfits, out[m1[srt]], repstr(outfile,'_all',''), clobber=clobber
    im_mwrfits, outseds[m1[srt]], repstr(outsedsfile,'_all',''), clobber=clobber

stop
return
end
    
