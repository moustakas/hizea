;pro build_witt, witt_dust, wittfile=wittfile
;; build the A(lambda) vs lambda look-up table using the Witt & Gordon
;; (2000) dust models and the flux-ratio method of Gordon et al. (2000)
;
;; choose three metallicities four ages, and five a_d values 
;    read_many_fit, dust_info, sf_type='c' ; continuous star formation
;    metal = [0.0,0.4]
;    ages = [10.0,50.0,100.0]
;;   metal = [-0.4,0.0,0.4]
;;   ages = [10.0,50.0,100.0,1000.0]
;    ad = [0.01,0.25,0.5,0.75,0.95]
;    ndust = n_elements(metal)*n_elements(ages)*n_elements(ad)
;
;    irx = range(0.1,1D6,1D3,/log)
;    nirx = n_elements(irx)
;    witt_dust = {irx: float(irx), a1500: float(irx*0.0), $
;      a1500_grid: fltarr(nirx,ndust), params: strarr(ndust)}
;
;    counter = 0
;    for im = 0, n_elements(metal)-1 do begin
;       for ia = 0, n_elements(ages)-1 do begin
;          for id = 0, n_elements(ad)-1 do begin
;             witt_dust.params[counter] = 'Z'+repstr(string(metal[im],format='(F4.1)'),' ','+')+$
;               '_Age'+string(ages[ia],format='(I4.4)')+'_ad'+string(ad[id],format='(F4.2)')
;             for jj = 0L, nirx-1L do witt_dust.a1500_grid[jj,counter] = fr2atten(irx[jj],$
;               '',1500.0,ages[ia],metal[im],ad[id],dust_info,/silent)
;             counter++
;          endfor
;       endfor
;    endfor
;
;; compute the media A(1500)-IRX relation and write out
;    witt_dust.a1500 = total(witt_dust.a1500_grid,2)/float(ndust)
;;   for jj = 0L, nirx-1L do witt_dust.a1500[jj] = djs_median(witt_dust.a1500_grid[jj,*])
;
;    im_mwrfits, witt_dust, wittfile, /clobber
;
;return
;end

pro sigmasfr_lir, clobber=clobber, rebuild_witt=rebuild_witt
; jm12mar19ucsd - compute L(IR) for the SIGMASFR sample

    path = sigmasfr_path()

    outfile = path+'sigmasfr_lir.fits'
    wittfile = path+'sigmasfr_witt.fits'
    
    kcorr = mrdfits(path+'sigmasfr_kcorrect.fits.gz',1)
    cat = mrdfits(path+'sigmasfr_photometry.fits.gz',1)
    ngal = n_elements(cat)

    out = struct_addtags(struct_trimtags(cat,select=['galaxy','z']),$
      replicate({$
      lir_chary: -999.0, lir_dale: -999.0, lir_rieke: -999.0, $
      sfr_chary: -999.0, sfr_dale: -999.0, sfr_rieke: -999.0, $
;     uvsfr_chary: -999.0, uvsfr_dale: -999.0, uvsfr_rieke: -999.0, $
      indx_chary: -999.0, indx_dale: -999.0, indx_rieke: -999.0, $
      irx_chary: -999.0, irx_dale: -999.0, irx_rieke: -999.0, $
      a1500_chary: -999.0, a1500_dale: -999.0, a1500_rieke: -999.0},ngal))
    
    filters = strtrim(sigmasfr_filterlist(),2)
    these = where(filters eq 'wise_w3.par' or filters eq 'wise_w4.par')

    good = where(total(cat.maggies[these] gt 0,1) gt 0)
    out[good].lir_chary = alog10(im_wise2lir(cat[good].z,cat[good].maggies[these],$
      cat[good].ivarmaggies[these],/chary,model_indx=indx_chary,debug=0))
stop
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

;    for i=0L, n_elements(good)-1L do begin
;      SFR = wr11_lir_to_sfr(wr11_f24_to_lir($
;        flux=cat[good[i]].maggies[10]*3631.d3, $
;        zobj=cat[good[i]].z))
;      splog, sfr, out[good[i]].sfr_chary
;    endfor

; build the Witt lookup table
    if file_test(wittfile+'.gz') eq 0 or keyword_set(rebuild_witt) then $
      build_witt, witt, wittfile=wittfile else $
        witt = mrdfits(wittfile+'.gz',1)

    ;l1500 = kcorr[good].k_uvflux[0]*(4.0*!dpi*dluminosity(out[good].z,/cm)^2)/im_lsun()
    l1500 = kcorr[good].k_uvflux[0]*(4.0*!dpi*dluminosity(out[good].z,/cm)^2)*1500./im_lsun()
    out[good].irx_rieke = 10.0^out[good].lir_rieke/l1500
    out[good].irx_dale = 10.0^out[good].lir_dale/l1500
    out[good].irx_chary = 10.0^out[good].lir_chary/l1500

    linterp, witt.irx, witt.a1500, out[good].irx_rieke, a1500_rieke, missing=0.0
    linterp, witt.irx, witt.a1500, out[good].irx_dale, a1500_dale, missing=0.0
    linterp, witt.irx, witt.a1500, out[good].irx_chary, a1500_chary, missing=0.0

    out[good].a1500_rieke = a1500_rieke
    out[good].a1500_dale = a1500_dale
    out[good].a1500_chary = a1500_chary

; get UV-based SFRs
;   out[good].uvsfr_chary = 1.4D-28*l1500*im_lsun()*1500/im_light(/ang)*10.0^(0.4*out[good].a1500_chary) 
    
; write out    
    im_mwrfits, out, outfile, clobber=clobber
stop
return
end
    
