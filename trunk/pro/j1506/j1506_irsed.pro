pro oplot_irsed, model, z=z, indx=indx, color=color, line=line, submm=submm
    nmodel = n_elements(model.lir)
    dlum = dluminosity(z,/cm)
    modelflam = interpolate(model.flux,indx,/grid)/(4.0*!dpi*dlum^2)/(1.0+z)
    good = where(modelflam gt 0.0,ngood)
    modelwave = model.wave[good]*(1.0+z)
;   modelfnu = modelflam[good]*rebin(reform(model.wave[good],ngood,1),ngood,nmodel)^2/im_light(/ang)
    modelfnu = modelflam[good]*rebin(reform(modelwave,ngood,1),ngood,nmodel)^2/im_light(/ang)*1D26 ; [cgs->mJy]

    modelmab = -2.5*alog10(modelfnu)+16.4
    modelwave = modelwave/1D4 ; [micron]
;   djs_oplot, modelwave, modelfnu, color=im_color(color), line=line, thick=6
    djs_oplot, modelwave, modelmab, color=im_color(color), line=line, thick=6

; print the 450 and 850 micron flux density
    submm = interpol(modelfnu,modelwave,[450.0,850.0])
return
end

pro j1506_irsed
; jm12feb27ucsd - analyze the IR SED of J1506

    j1506path = j1506_path()
    cat = mrdfits(sigmasfr_path()+'sigmasfr_photometry.fits.gz',1)
    cat = cat[where(strmatch(cat.galaxy,'*J1506+54*'))]

    filters = strtrim(sigmasfr_filterlist(),2)
    iwise = where(strmatch(filters,'*wise*',/fold))
    weff = k_lambda_eff(filterlist=wise_filterlist())/1D4
    
    lir_chary = im_wise2lir(cat.z,cat.maggies[iwise[2:3]],cat.ivarmaggies[iwise[2:3]],$
      /chary,model_indx=indx_chary,chi2=chi2_chary) ;,/debug)
    lir_rieke = im_wise2lir(cat.z,cat.maggies[iwise[2:3]],cat.ivarmaggies[iwise[2:3]],$
      /rieke,model_indx=indx_rieke,chi2=chi2_rieke);,/debug)
    lir_dale = im_wise2lir(cat.z,cat.maggies[iwise[2:3]],cat.ivarmaggies[iwise[2:3]],$
      /dale,model_indx=indx_dale,chi2=chi2_dale);/debug)
    lir = [lir_chary,lir_rieke,lir_dale]
    sfr = lir*4.5D-44*im_lsun()

; make the QAplot    
    chary = read_01chary()
    rieke = read_09rieke()
    dale = read_02dale()
    
    psfile = j1506path+'j1506_irsed.ps'
    yrange = [20.0,8.0]
    im_plotconfig, 0, pos, psfile=psfile, height=4.5, width=5.6, xmargin=[1.2,1.7]
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=9, $
      xrange=[1,2000], yrange=yrange, /xlog, $
      xtitle='Observed-Frame Wavelength (\mu'+'m)', $
      ytitle='AB magnitude'
    axis, /yaxis, ytitle='Flux Density (mJy)', yrange=10.0^(-0.4*(yrange-16.4)), $
      ysty=1, /ylog
;   djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;     xrange=[1,1000], yrange=[1D-2,300], /xlog, /ylog, $
;     xtitle='Observed-Frame Wavelength (\mu'+'m)', $
;     ytitle='Flux Density (mJy)'
    im_legend, ['Chary & Elbaz','Dale & Helou','Rieke'], /left, /top, box=0, $
      charsize=1.6, color=['orange','dodger blue','forest green'], $
      line=[0,5,3], thick=8

    oplot_irsed, chary, indx=indx_chary, z=cat.z, color='orange', line=0, submm=chary_submm
    oplot_irsed, dale, indx=indx_dale, z=cat.z, color='dodger blue', line=5, submm=dale_submm
    oplot_irsed, rieke, indx=indx_rieke, z=cat.z, color='forest green', line=3, submm=rieke_submm

    djs_oplot, 450*[1,1], [15,10], thick=5
    djs_oplot, 850*[1,1], [16,11], thick=5
    
    fmJy = cat.maggies[iwise]*10^(+0.4*16.4)
    emJy = 10^(+0.4*16.4)/sqrt(cat.ivarmaggies[iwise])
    mab = maggies2mag(cat.maggies[iwise],ivarmaggies=cat.ivarmaggies[iwise],magerr=maberr)
;   niceprint, weff, mab, maberr
    
;   oploterror, weff, fmJy, emJy, psym=symcat(16), $
;     symsize=2.0, errthick=6
    oploterror, weff, mab, maberr, psym=symcat(16), $
      symsize=2.0, errthick=6

    im_plotconfig, psfile=psfile, /psclose, /pdf

    out = replicate({model: '', chi2: 0.0, lir: 0.0, sfr: 0.0, f450_mJy: 0.0, f850_mJy: 0.0},3)
    out.model = ['chary','rieke','dale']
    out.lir = lir
    out.chi2 = [chi2_chary,chi2_rieke,chi2_dale]
    out.sfr = sfr
    out.f450_mJy = [chary_submm[0],rieke_submm[0],dale_submm[0]]
    out.f850_mJy = [chary_submm[1],rieke_submm[1],dale_submm[1]]
    struct_print, out
    
    
stop
    
return
end
    
