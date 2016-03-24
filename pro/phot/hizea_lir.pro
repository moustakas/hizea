pro hizea_lir, debug=debug, clobber=clobber
; jm15nov17siena - compute L(IR) for the full Sample 2 sample

    photdir = hizea_path(/phot)

    outfile = photdir+'hizea_sample2_lir.fits'
    outsedsfile = photdir+'hizea_sample2_lir_seds.fits'
    
;   kcorr = mrdfits(path+'massprofiles_kcorrect.fits.gz',1)
    cat = mrdfits(photdir+'hizea_sample2_phot.fits.gz',1)
;   cat = cat[0:10]
    ngal = n_elements(cat)

    gal = strmid(cat.galaxy,0,5)+strmid(cat.galaxy,10,1)+$
      strmid(cat.galaxy,12,4)
    
    light = 2.99792458D18       ; speed of light [A/s]
    filters = strtrim(hizea_filterlist(),2)
    keep = where(strmatch(filters,'*wise*') or strmatch(filters,'*irac*'))
    
    filters = filters[keep]
    these = where(filters eq 'wise_w3.par' or filters eq 'wise_w4.par')
    weff = k_lambda_eff(filterlist=filters)

; nuLnu is the rest-frame luminosity 
    out = struct_addtags(struct_trimtags(cat,select=['galaxy','z']),$
      replicate({$
      lir_chary: -999.0, lir_dale: -999.0, $
;     lir_rieke: -999.0, $
      nuLnu_chary: dblarr(3)-999.0, nuLnu_dale: dblarr(3)-999.0, $
;     nuLnu_rieke: dblarr(3)-999.0, $
      ivar_nulnu_chary: dblarr(3)-999.0, ivar_nulnu_dale: dblarr(3)-999.0, $
;     ivar_nulnu_rieke: dblarr(3)-999.0, $
      sfr_chary: -999.0, sfr_dale: -999.0, $
;     sfr_rieke: -999.0, $
;     uvsfr_chary: -999.0, uvsfr_dale: -999.0, uvsfr_rieke: -999.0, $
      indx_chary: -999.0, indx_dale: -999.0},ngal))
;     indx_rieke: -999.0, $
;     irx_chary: -999.0, irx_dale: -999.0, irx_rieke: -999.0, $
;     a1500_chary: -999.0, a1500_dale: -999.0, a1500_rieke: -999.0},ngal))

    outseds = replicate({$
      modelwave_chary: fltarr(1366), modelwave_dale: fltarr(1930), modelwave_rieke: fltarr(6270), $
      modelmab_chary: fltarr(1366), modelmab_dale: fltarr(1930), modelmab_rieke: fltarr(6270)},$
      ngal)
    
; throw out upper limits to make sure the K-corrections in
; im_simple_kcorrect get computed from the band where the object was
; detected 
    good = where(total(cat.maggies[keep[these]] gt 0,1) gt 0)
    wmaggies = cat[good].maggies[keep[these]]
;   wivarmaggies = cat[good].ivarmaggies[keep[these]]
    wivarmaggies = cat[good].ivarmaggies[keep[these]]*(wmaggies gt 0)
    
    out[good].lir_chary = alog10(im_wise2lir(cat[good].z,wmaggies,wivarmaggies,$
      /chary,model_indx=indx_chary,nulnu=nulnu_chary,debug=debug,$
      ivar_nulnu=ivar_nulnu_chary,modelwave=modelwave_chary,modelmab=modelmab_chary))
    out[good].lir_dale = alog10(im_wise2lir(cat[good].z,wmaggies,wivarmaggies,$
      /dale,model_indx=indx_dale,nulnu=nulnu_dale,$
      ivar_nulnu=ivar_nulnu_dale,modelwave=modelwave_dale,modelmab=modelmab_dale))
;   out[good].lir_rieke = alog10(im_wise2lir(cat[good].z,wmaggies,wivarmaggies,$
;     /rieke,model_indx=indx_rieke,nulnu=nulnu_rieke,$
;     ivar_nulnu=ivar_nulnu_rieke,modelwave=modelwave_rieke,modelmab=modelmab_rieke))

    outseds[good].modelwave_chary = modelwave_chary
    outseds[good].modelmab_chary = modelmab_chary
    outseds[good].modelwave_dale = modelwave_dale
    outseds[good].modelmab_dale = modelmab_dale
;   outseds[good].modelwave_rieke = modelwave_rieke
;   outseds[good].modelmab_rieke = modelmab_rieke

    out[good].indx_chary = indx_chary
    out[good].indx_dale = indx_dale
;   out[good].indx_rieke = indx_rieke

    out[good].nulnu_chary = nulnu_chary
    out[good].nulnu_dale = nulnu_dale
;   out[good].nulnu_rieke = nulnu_rieke

    out[good].ivar_nulnu_chary = ivar_nulnu_chary
    out[good].ivar_nulnu_dale = ivar_nulnu_dale
;   out[good].ivar_nulnu_rieke = ivar_nulnu_rieke

    fact = 4.5D-44*im_lsun()
    out[good].sfr_chary = 10.0^out[good].lir_chary*fact
    out[good].sfr_dale = 10.0^out[good].lir_dale*fact
;   out[good].sfr_rieke = 10.0^out[good].lir_rieke*fact

;    for i=0L, n_elements(good)-1L do begin
;      SFR = wr11_lir_to_sfr(wr11_f24_to_lir($
;        flux=cat[good[i]].maggies[10]*3631.d3, $
;        zobj=cat[good[i]].z))
;      splog, sfr, out[good[i]].sfr_chary
;    endfor

; write out    
    im_mwrfits, out, outfile, clobber=clobber
    im_mwrfits, outseds, outsedsfile, clobber=clobber

; build a QAplot
    psfile = photdir+'qa_hizea_sample2_lir.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, $
      xmargin=[1.8,0.3], width=6.4, thick=8
    for ii = 0, ngal-1 do begin
       wave_chary = outseds[ii].modelwave_chary
       fnu_chary = 10^(-0.4*(outseds[ii].modelmab_chary+48.6))
       maggies_chary = reform(k_project_filters(k_lambda_to_edges(wave_chary),$
         fnu_chary*light/wave_chary^2,filterlist=filters))*10^(-0.4*48.6)*1D26 ; [mJy]

       wave_dale = outseds[ii].modelwave_dale
       fnu_dale = 10^(-0.4*(outseds[ii].modelmab_dale+48.6))
       maggies_dale = reform(k_project_filters(k_lambda_to_edges(wave_dale),$
         fnu_dale*light/wave_dale^2,filterlist=filters))*10^(-0.4*48.6)*1D26 ; [mJy]

       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         /xlog, /ylog, xrange=[1.0,5E3], yrange=[0.001,3000], $
         xtitle='Wavelength (\mu'+'m)', ytitle='F_{\nu} (mJy)'
       djs_oplot, wave_chary/1D4, fnu_chary*1D26
       djs_oplot, wave_dale/1D4, fnu_dale*1D26, line=5

       djs_oplot, weff/1D4, maggies_chary, psym=6, symsize=3, color='orange'
       djs_oplot, weff/1D4, maggies_dale, psym=cgsymcat(9), $
         symsize=3, color=cgcolor('forest green')
       djs_oplot, weff/1D4, cat[ii].maggies[keep]*10^(-0.4*48.6)*1D26, $
         psym=7, symsize=3, color='blue'

       im_legend, [gal[ii],'z='+strmid(strtrim(cat[ii].z,2),0,5)], $
         /left, /top, box=0
       im_legend, ['Chary & Elbaz (2001)','Dale & Helou (2002)'], /right, /bottom, $
         box=0, line=[0,5], pspacing=1.9, charsize=1.5
       im_legend, ['WISE observed','CE01 model','DH02 model'], /right, /top, $
         box=0, charsize=1.5, psym=[7,6,9], color=['blue','orange','forest green']

;      if strmatch(gal,'*1613') then plots, 33D9
;      57.5microJy at 33 GHz
       
;      cc = get_kbrd(1)
    endfor
    im_plotconfig, psfile=psfile, /psclose, /pdf


stop
    
    
return
end
