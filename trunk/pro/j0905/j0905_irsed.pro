pro oplot_irsed, model, z=z, indx=indx, color=color, line=line
    nmodel = n_elements(model.lir)
    dlum = dluminosity(z,/cm)
    modelflam = interpolate(model.flux,indx,/grid)/(4.0*!dpi*dlum^2)/(1.0+z)
    good = where(modelflam gt 0.0,ngood)
    modelwave = model.wave[good]*(1.0+z)
    modelmab = -2.5*alog10(modelflam[good]*rebin(reform(model.wave[good],ngood,1),ngood,nmodel)^2/im_light(/ang))-48.6
    djs_oplot, modelwave/1D4, modelmab, color=im_color(color), line=line, thick=6
return
end

pro j0905_irsed
; jm12feb27ucsd - analyze the IR SED of J0905

    j0905path = j0905_path()
    cat = mrdfits(j0905path+'j0905_photometry.fits.gz',1)

    mm = cat.maggies[9:12]
    ii = cat.ivarmaggies[9:12]
    weff = k_lambda_eff(filterlist=wise_filterlist())
    
    lir_chary = im_wise2lir(cat.z,mm[2:3],ii[2:3],/chary,model_indx=indx_chary)
    lir_rieke = im_wise2lir(cat.z,mm[2:3],ii[2:3],/rieke,model_indx=indx_rieke)
    lir_dale = im_wise2lir(cat.z,mm[2:3],ii[2:3],/dale,model_indx=indx_dale)
    lir = [lir_chary,lir_rieke,lir_dale]
    sfr = lir*4.5D-44*im_lsun()
    niceprint, ['chary','rieke','dale'], lir, sfr

;; get the total attenuation
;    kk = mrdfits(j0905path+'j0905_kcorrect.fits.gz',1)
;    l1500 = kk.k_uvflux[0]*(4.0*!dpi*dluminosity(cat.z,/cm)^2)
;    irx = lir/l1500*im_lsun()
;    linterp, witt.irx, witt.a1500, irx, a1500, missing=0.0

; make the QAplot    
    chary = read_01chary()
    rieke = read_09rieke()
    dale = read_02dale()
    
    psfile = j0905path+'j0905_irsed.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.5
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xrange=[1,500], yrange=[24,10], /xlog, $
      xtitle='Observed-Frame Wavelength (\mu'+'m)', $
      ytitle='AB Magnitude'
    im_legend, ['Chary & Elbaz','Dale & Helou','Rieke'], /left, /top, box=0, $
      charsize=1.6, color=['orange','dodger blue','forest green'], $
      line=[0,5,3], thick=8

    oplot_irsed, chary, indx=indx_chary, z=cat.z, color='orange', line=0
    oplot_irsed, dale, indx=indx_dale, z=cat.z, color='dodger blue', line=5
    oplot_irsed, rieke, indx=indx_rieke, z=cat.z, color='forest green', line=3

    mab = maggies2mag(mm,ivarmaggies=ii,magerr=maberr)
    oploterror, weff/1D4*(1+cat.z), mab, maberr, psym=symcat(16), $
      symsize=2.0, errthick=6

    im_plotconfig, psfile=psfile, /psclose, /pdf, /pskeep

return
end
    
