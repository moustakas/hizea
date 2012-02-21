pro j0905_render_postplot, xx, pos, xrange=xrange, yrange=yrange, $
  xtitle=xtitle, ytitle=ytitle, binsize=binsize, noerase=noerase, $
  xtickinterval=xtickinterval, nomedian=nomedian, monte=monte

    yrange = [0,1.4]
;   yrange = [0,1.05]
    if (n_elements(xrange) eq 0) then xrange = minmax(xx)*[0.9,1.1]
    if (n_elements(binsize) eq 0) then $
      binsize = (xrange[1]-xrange[0])/ceil(0.3*sqrt(n_elements(xx)))

    plot, [0], [0], xsty=5, ysty=5, /nodata, position=pos, $
      yrange=yrange, xrange=xrange, noerase=noerase
    im_plothist, xx, bin=binsize, /peak, /overplot, /fill, $
      fcolor=im_color('grey80'), xhist, yhist
;   if keyword_set(nomedian) eq 0 then $
;     djs_oplot, median(xx)*[1,1], !y.crange, line=5, thick=8
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, position=pos, $
      yrange=yrange, xrange=xrange, ytitle=ytitle, $
      xtitle=textoidl(xtitle), xtickinterval=xtickinterval, $
      ytickname=replicate(' ',10)
    if n_elements(monte) ne 0 then begin
       im_plothist, monte, bin=binsize*1.0, mxhist, myhist, /noplot
       im_plothist, monte, bin=binsize*1.0, /overplot, line=1, $
         normfactor=max(myhist)/(yrange[1]*0.95)
;        normfactor=max(myhist)/max(yhist)/1.2
    endif

return
end

pro j0905_plots, supergrid, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber
; jm12feb14ucsd - build plots for the paper

    common j0905_post, dist
    
    j0905path = j0905_path()
    isedpath = j0905_path(/isedfit)
    isedfit_sfhgrid_dir = j0905_path(/monte)
    sfhgrid_paramfile = getenv('HIZEA_DIR')+'/pro/j0905/j0905_sfhgrid.par'

    super = get_j0905_supergrid(supergrid,nsuper=nsuper)
    struct_print, super

    filters = j0905_filterlist()
    weff = k_lambda_eff(filterlist=filters)
    toflam = 10^(-0.4*48.6);*im_light(/ang)/weff^2
;   toflam = 10^(-0.4*48.6)*im_light(/ang)/weff^2
    ndraw = isedfit_ndraw()

; gather the photometry and restore the results
    cat = mrdfits(j0905path+'j0905_photometry.fits.gz',1)
;   cat = mrdfits(j0905path+'j0905+5759_isedfit_input_v2.fits.gz',1)

    if (n_elements(dist) eq 0) then begin
       paramfile = isedpath+'j0905_supergrid05_isedfit.par' ; smc dust
;      paramfile = isedpath+'j0905_supergrid04_isedfit.par' ; smc dust
       mstar = j0905_reconstruct_posterior(paramfile,post=post,$
         isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,iopath=isedpath,$
         age=age,Z=Z,tau=tau,sfr0=sfr0,av=av,sfrage=sfrage,$
         chunkindx=chunkindx,modelindx=modelindx,indxage=ageindx,$
         trunctau=tautrunc,tburst=tburst,fburst=fburst,dtburst=dtburst,$
         sfrpeak=sfrpeak,timesinceburst=timesinceburst)

       dist = {mstar: fltarr(ndraw), tau: fltarr(ndraw), age: fltarr(ndraw), $
         av: fltarr(ndraw), Z: fltarr(ndraw), tautrunc: fltarr(ndraw), $
         tburst: fltarr(ndraw), fburst: fltarr(ndraw), dtburst: fltarr(ndraw), $
         sfrage: fltarr(ndraw), sfr0: fltarr(ndraw), sfrpeak: fltarr(ndraw), $
         timesinceburst: fltarr(ndraw)}

       dist.mstar = mstar
       dist.tau = tau
       dist.age = age
       dist.av = av
       dist.Z = Z
       dist.tautrunc = tautrunc
       dist.tburst = tburst
       dist.fburst = fburst
       dist.dtburst = dtburst
       dist.sfrage = sfrage
       dist.sfr0 = sfr0
       dist.sfrpeak = sfrpeak
       dist.timesinceburst = timesinceburst
    endif

;; reconstruct the posterior models       
;       temp = replicate({zobj: cat.z, chi2: 0.0, chunkindx: 0L, $
;         modelindx: 0L, ageindx: 0L, scale: 0.0},ndraw)
;       temp.chunkindx = chunkindx
;       temp.modelindx = modelindx
;       temp.ageindx = ageindx
;       temp.scale = post.scale
;       postmodel = isedfit_restore(paramfile,in_isedfit=temp,$
;         isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,iopath=isedpath)
;;      plot, postmodel[0].wave, postmodel[0].flux, /xlog, xr=[1E4,4E4], yr=[30,20], ps=3
;;      for ii = 1, ndraw-1 do djs_oplot, postmodel[ii].wave, postmodel[ii].flux, ps=3
;    endif

; ---------------------------------------------------------------------------
; posterior distributions plots
    monte = mrdfits(isedfit_sfhgrid_dir+'sfhgrid02/fsps/smc/chab_montegrid.fits.gz',1)
;   monte = mrdfits(isedfit_sfhgrid_dir+'sfhgrid01/bc03/smc/chab_montegrid.fits.gz',1)

; mass, reddening, and metallicity
    psfile = j0905path+'j0905_massdust_posteriors.eps'
    im_plotconfig, 5, pos, psfile=psfile, xspace=0.3, yspace=1.0, $
      xmargin=[1.0,0.5], height=3.0*[1,1], width=3.5*[1,1], charsize=1.8

    j0905_render_postplot, dist.mstar, pos[*,0], xrange=[10,11.3], $
      ytitle='Probability', xtitle='log (M_{*}/M'+sunsymbol()+')', $
      xtickinterval=0.5

    j0905_render_postplot, dist.Z/0.019, pos[*,1], xrange=[0.5,1.7], $
      xtitle='Z/Z'+sunsymbol(), /noerase, /nomedian, monte=monte.z/0.019, $
      xtickinterval=0.5
    
    j0905_render_postplot, dist.av, pos[*,2], xrange=minmax(monte.av)+[-0.05,0.05], $
      ytitle='Probability', xtitle='A_{V} (mag)', /noerase, monte=monte.av
    
    im_plotconfig, /psclose, psfile=psfile, /pdf, /pskeep

; SFH posteriors
    psfile = j0905path+'j0905_sfh_posteriors.eps'
    im_plotconfig, 17, pos, psfile=psfile, xspace=[0.3,0.3], yspace=1.0, $
      xmargin=[1.0,0.5], height=2.1*[1,1], width=2.7*[1,1,1], charsize=1.4

    j0905_render_postplot, dist.sfr0, pos[*,0], $
      xrange=[1.8,3.0], $      ; xrange=minmax(dist.sfr0-logmu)*[0.99,1.01], $
;     xrange=[-0.3,2.2], $      ; xrange=minmax(dist.sfr0-logmu)*[0.99,1.01], $
      ytitle='Probability', xtitle='log (SFR /M'+sunsymbol()+' yr^{-1})' ;, monte=dist.bigsfr+10
    
    j0905_render_postplot, dist.sfrpeak, pos[*,1], /noerase, $
      xtitle='log (SFR_{peak} /M'+sunsymbol()+' yr^{-1})'
    
    j0905_render_postplot, dist.tau, pos[*,2], /noerase, xrange=[-0.2,7.2], $
      xtitle='\tau (Gyr)'

    j0905_render_postplot, dist.tautrunc*1E3, pos[*,3], /noerase, $
      ytitle='Probability', xtitle='\tau_{trunc} (Myr)', xrange=[0,55]
    
    j0905_render_postplot, dist.dtburst*1E3, pos[*,4], /noerase, $
      xtitle='\sigma_{burst} (Myr)', xrange=[0,105]

    j0905_render_postplot, dist.timesinceburst*1E3, pos[*,5], /noerase, $
      xtitle='\Delta'+'t (Myr)', xrange=[-249,130], xtickinterval=100

    im_plotconfig, /psclose, psfile=psfile, /pdf, /pskeep

; ---------------------------------------------------------------------------
; FSPS SED plot for the paper
    paramfile = isedpath+'j0905_supergrid05_isedfit.par'
    model = isedfit_restore(paramfile,ised,iopath=isedpath,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir);,/fnu)

    psfile = j0905path+'j0905_fsps_smc.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=4.8, ymargin=[1.0,1.0], $
      xmargin=[1.4,0.4], width=6.7, charthick=3.0

    xtitle1 = textoidl('Observed-Frame Wavelength (\AA)')
    xtitle2 = textoidl('Rest-Frame Wavelength (\AA)')
    ytitle1 = textoidl('Apparent AB Magnitude')

    yrange = [25.0,17.5]
    xrange1 = [1000,80000]
;   xrange1 = [1000,60000]

    ticks1 = [2000,4000,12000,30000,60000]
    ticks2 = [1000,2000,4000,8000,16000,70000]

    xrange2 = xrange1/(1.0+ised.zobj)
    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
      xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, xlog=1, $ ; /ylog, 
      position=pos, xtickformat='(I0)', xticks=n_elements(ticks1)-1, xtickv=ticks1
    axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2, xlog=1, $
      xticks=n_elements(ticks2)-1, xtickv=ticks2, xtickformat='(I0)'
    djs_oplot, model.wave, model.flux, line=0, color='black'
    
    used = where((ised[0].maggies gt 0.0) and $ ; used in the fitting
      (ised[0].ivarmaggies gt 0.0),nused)
    notused = where((ised.maggies gt 0.0) and $
      (ised.ivarmaggies eq 0.0),nnotused)

    mab = -2.5*alog10(ised.bestmaggies[used])
;   djs_oplot, weff[used], mab, psym=symcat(6,thick=6), $
;     symsize=2.5, color=im_color('steel blue')

    if (nused ne 0L) then begin
       mab = maggies2mag(ised.maggies[used],$
         ivar=ised[0].ivarmaggies[used],magerr=mab_err)
       oploterror, weff[used], mab, mab_err, psym=symcat(16), $
         symsize=1.5, color=im_color('firebrick'), $
         errcolor=im_color('firebrick'), errthick=8
    endif

    if (nnotused ne 0) then begin
       mab = maggies2mag(ised.maggies[notused],$
         ivar=ised.ivarmaggies[notused],magerr=mab_err)
       oploterror, weff[notused], mab, mab_err, psym=symcat(16), $
         symsize=2.0, color=im_color('forest green'), $
         errcolor=im_color('forest green'), errthick=8
    endif
   
    im_plotconfig, /psclose, psfile=psfile, /pdf, /pskeep

stop    
    
    
    
; ---------------------------------------------------------------------------
; plot the maximum likelihood sfh
    ised = mrdfits(isedpath+'j0905_bc03_chab_smc_sfhgrid01.fits.gz',1)
    
    psfile = j0905path+'j0905_bestfit_sfh.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0

    sfh = alog10(isedfit_reconstruct_sfh(ised,outage=outage,dage=0.001D,nagemax=10000))+alog10(ised.scale)
    
    yrange = [-1,3]
    xrange = [5,ised.age]
    xtitle = 'Time (Gyr)'
    ytitle = textoidl('log (SFR / M'+sunsymbol()+' yr^{-1})')

    plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
      xsty=1, ysty=3, xtitle=xtitle, ytitle=ytitle, position=pos
    djs_oplot, outage, sfh, line=0, thick=8;, color='pink', psym=6

    im_legend, [$
      '\tau = '+strtrim(string(ised.tau,format='(F12.1)'),2)+' Gyr', $
      't_{burst} = '+strtrim(string(ised.tburst,format='(F12.2)'),2)+' Gyr',$
      '\sigma_{burst} = '+strtrim(string(ised.dtburst*1E3,format='(F12.1)'),2)+' Myr',$
      '\tau_{trunc} = '+strtrim(string(ised.tautrunc*1E3,format='(F12.1)'),2)+' Myr'], $
      /left, /top, box=0, charsize=1.6
    
    im_plotconfig, /psclose, psfile=psfile, /pdf, /pskeep

stop    
    
; ---------------------------------------------------------------------------
; compare the LRIS spectrum to the best-fit model
    paramfile = isedpath+'j0905_supergrid04_isedfit.par'
    model = isedfit_restore(paramfile,ised,iopath=isedpath,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,/flam)
    lris = mrdfits(j0905path+'spectra/j0905+5759_lris_flux_v120216.fits',1)
    
    psfile = j0905path+'j0905_lriscompare.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=4.8, ymargin=[1.0,1.0], $
      xmargin=[1.4,0.4], width=6.7, charthick=3.0

    xtitle1 = textoidl('Observed-Frame Wavelength (\AA)')
    xtitle2 = textoidl('Rest-Frame Wavelength (\AA)')
    ytitle1 = textoidl('Flux (10^{-17} '+flam_units()+')')

    scale = 1D17
    yrange = [0.05,15]
    xrange1 = [3000,11000]

    xrange2 = xrange1/(1.0+ised.zobj)
    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
      xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, position=pos
    axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2
    im_legend, ['LRIS','BC03 Model'], /right, /top, box=0, $
      charsize=1.6, line=0, pspacing=1.9, color=['grey','black']
    
    djs_oplot, model.wave, scale*model.flux, line=0, color='black'
    djs_oplot, lris.wavelength, smooth(scale*lris.flux,5), color='grey'
    
    im_plotconfig, /psclose, psfile=psfile, /pdf, /pskeep

stop    
    
; ---------------------------------------------------------------------------
; compare the best-fitting SEDs
    psfile = j0905path+'j0905_dustcompare.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=4.8, ymargin=[1.0,1.0], $
      xmargin=[1.4,0.4], width=6.7, charthick=3.0

    xtitle1 = textoidl('Observed-Frame Wavelength (\AA)')
    xtitle2 = textoidl('Rest-Frame Wavelength (\AA)')
    ytitle1 = textoidl('Apparent AB Magnitude')

    yrange = [25,17.5]
    xrange1 = [1000,60000]

    ticks1 = [2000,4000,12000,30000]
    ticks2 = [1000,2000,4000,8000,16000]

    xrange2 = xrange1/1.7116
    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
      xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, xlog=1, $ ; /ylog, 
      position=pos, xtickformat='(I0)', xticks=n_elements(ticks1)-1, xtickv=ticks1
    axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2, xlog=1, $
      xticks=n_elements(ticks2)-1, xtickv=ticks2, xtickformat='(I0)'

    color = ['forest green','dodger blue','orange','black']
    for ii = 0, 3 do begin
       paramfile = isedpath+'j0905_supergrid0'+strtrim(ii+1,2)+'_isedfit.par'
       model = isedfit_restore(paramfile,ised,iopath=isedpath,$
         isedfit_sfhgrid_dir=isedfit_sfhgrid_dir)
       djs_oplot, model.wave, model.flux, line=0, color=im_color(color[ii])
    endfor
    im_legend, ['Calzetti','Charlot & Fall','Milky Way','SMC'], $
      /right, /bottom, box=0, line=0, color=color, pspacing=1.9, $
      charsize=1.7
    
    used = where((ised[0].maggies gt 0.0) and $ ; used in the fitting
      (ised[0].ivarmaggies gt 0.0),nused)
    notused = where((ised.maggies gt 0.0) and $
      (ised.ivarmaggies eq 0.0),nnotused)

;   mab = -2.5*alog10(ised.bestmaggies[used])
;   djs_oplot, weff[used], mab, psym=symcat(6,thick=6), $
;     symsize=2.5, color=im_color('steel blue')

    if (nused ne 0L) then begin
       mab = maggies2mag(ised.maggies[used],$
         ivar=ised[0].ivarmaggies[used],magerr=mab_err)
       oploterror, weff[used], mab, mab_err, psym=symcat(16), $
         symsize=1.5, color=im_color('firebrick'), $
         errcolor=im_color('firebrick'), errthick=8
    endif

    if (nnotused ne 0) then begin
       mab = maggies2mag(ised.maggies[notused],$
         ivar=ised.ivarmaggies[notused],magerr=mab_err)
       djs_oplot, weff[notused], mab, mab_err, psym=symcat(16), $
         symsize=2.0, color=im_color('forest green'), $
         errcolor=im_color('forest green'), errthick=8
    endif
   
    im_plotconfig, /psclose, psfile=psfile, /pdf, /pskeep

; ---------------------------------------------------------------------------
; SED plot for the paper
    paramfile = isedpath+'j0905_supergrid04_isedfit.par'
    model = isedfit_restore(paramfile,ised,iopath=isedpath,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir);,/fnu)
    lris = mrdfits(j0905path+'spectra/j0905+5759_lris_flux_v120216.fits',1)
    good = where(lris.flux gt 0)
    lriswave = lris[good].wavelength
    lrisflux = -2.5*alog10(lris[good].flux*lriswave^2/im_light(/ang))-48.6
;   lris.flux = lris.flux*lris.wavelength^2/im_light(/ang)
    
    psfile = j0905path+'j0905_sed_smc.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=4.8, ymargin=[1.0,1.0], $
      xmargin=[1.4,0.4], width=6.7, charthick=3.0

    xtitle1 = textoidl('Observed-Frame Wavelength (\AA)')
    xtitle2 = textoidl('Rest-Frame Wavelength (\AA)')
    ytitle1 = textoidl('Apparent AB Magnitude')
;   ytitle1 = textoidl('Flux (10^{-28} '+fnu_units()+')')
;   ytitle1 = textoidl('Flux (10^{-17} '+flam_units()+')')

;   scale = 1D28
;   scale = 1D17
;   yrange = [0.05,50]
    yrange = [25,17.5]
    xrange1 = [1000,60000]

    ticks1 = [2000,4000,12000,30000]
    ticks2 = [1000,2000,4000,8000,16000]

    xrange2 = xrange1/(1.0+ised.zobj)
    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
      xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, xlog=1, $ ; /ylog, 
      position=pos, xtickformat='(I0)', xticks=n_elements(ticks1)-1, xtickv=ticks1
    axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2, xlog=1, $
      xticks=n_elements(ticks2)-1, xtickv=ticks2, xtickformat='(I0)'

    djs_oplot, model.wave, model.flux, line=0
;   djs_oplot, model.wave, scale*model.flux, line=0
    djs_oplot, lriswave, lrisflux, color='orange'
;   djs_oplot, lris.wavelength, scale*lris.flux, color='orange'
    
    used = where((ised[0].maggies gt 0.0) and $ ; used in the fitting
      (ised[0].ivarmaggies gt 0.0),nused)
    notused = where((ised.maggies gt 0.0) and $
      (ised.ivarmaggies eq 0.0),nnotused)

    mab = -2.5*alog10(ised.bestmaggies[used])
    djs_oplot, weff[used], mab, psym=symcat(6,thick=6), $
      symsize=2.5, color=im_color('steel blue')
;   best_flam = scale*ised.bestmaggies[used]*toflam[used]
;   djs_oplot, weff[used], best_flam, psym=symcat(6,thick=6), $
;     symsize=2.5, color=im_color('steel blue')

    if (nused ne 0L) then begin
       mab = maggies2mag(ised.maggies[used],$
         ivar=ised[0].ivarmaggies[used],magerr=mab_err)
       oploterror, weff[used], mab, mab_err, psym=symcat(16), $
         symsize=2.0, color=im_color('firebrick'), $
         errcolor=im_color('firebrick'), errthick=8
    endif

;   phot_flam = scale*ised.maggies[used]*toflam[used]
;   photerr_flam = scale*toflam[used]/sqrt(ised.ivarmaggies[used])
;   oploterror, weff[used], phot_flam, photerr_flam, psym=symcat(16), $
;     symsize=1.0, color=im_color('firebrick'), $
;     errcolor=im_color('firebrick'), errthick=8

    if (nnotused ne 0) then begin
       mab = maggies2mag(ised.maggies[notused],$
         ivar=ised.ivarmaggies[notused],magerr=mab_err)
       djs_oplot, weff[notused], mab, mab_err, psym=symcat(16), $
         symsize=2.0, color=im_color('forest green'), $
         errcolor=im_color('forest green'), errthick=8
;      phot_flam = scale*ised.maggies[notused]*toflam[notused]
;      djs_oplot, weff[notused], phot_flam, psym=symcat(15), $
;        symsize=2.0, color=im_color('forest green')
    endif
   
    im_plotconfig, /psclose, psfile=psfile, /pdf, /pskeep

stop    
    
; ---------------------------------------------------------------------------
; plot the NxN distributions of parameters as an upper triangle
    this = 0 ; z=9.56 solution
    
    psfile = j0905path+'j0905_manyd.eps'
    splog, 'Writing '+psfile
    manyd = transpose([$
      [dist.mstar],$
      [dist.Z/0.019],$
      [dist.av],$
      [dist.timesinceburst]*1E3,$
      [dist.tautrunc]*1E3,$
      [dist.sfr0],$
      [dist.sfrpeak]])
    label = textoidl([$
      'log (M/M'+sunsymbol()+')',$
      'Z/Z'+sunsymbol(),$
      'A_{V} (mag)',$
      'time (Myr)',$
      '\tau_{trunc} (Myr)',$
      'log (\psi_{max}/M'+sunsymbol()+' yr^{-1})',$
      'log (\psi_{0}/M'+sunsymbol()+' yr^{-1})'])
    
    im_manyd_scatterplot, fltarr(ndraw)+1, manyd, psfile, label=label, $
      axis_char_scale=1.4, /internal, outliers=1, $
      /nogrey, levels=errorf((dindgen(2)+1)/sqrt(2)), /upper
    spawn, 'ps2pdf '+psfile+' '+repstr(psfile,'.eps','.pdf'), /sh


; ---------------------------------------------------------------------------
; plot the NxN distributions of parameters as an upper triangle
    this = 0 ; z=9.56 solution
    
    psfile = datapath+'santorini_manyd.eps'
    splog, 'Writing '+psfile
    manyd = transpose([[mstar[*,this]],[sfrage[*,this]],[Z[*,this]/0.019],[av[*,this]],[sfr0[*,this]]])
    label = textoidl(['log (M/M'+sunsymbol()+')','t_{w} (Gyr)','Z/Z'+sunsymbol(),'A_{V} (mag)','log (\psi/M'+sunsymbol()+' yr^{-1})'])
    
    im_manyd_scatterplot, fltarr(ndraw)+1, manyd, psfile, label=label, $
      axis_char_scale=1.4, /internal, outliers=1, $
      /nogrey, levels=errorf((dindgen(2)+1)/sqrt(2)), /upper
    spawn, 'ps2pdf '+psfile, /sh

stop    

; ---------------------------------------------------------------------------
; posterior distributions: mass, sfrage, sfr100, av, Z, sSFR
    psfile = datapath+'santorini_posteriors.eps'
    im_plotconfig, 14, pos, psfile=psfile, height=5.0, ymargin=[1.0,1.0]

    plot, [0], [0], xsty=5, ysty=5, /noerase, /nodata, position=pos2[*,0], yrange=[0,1.05], $
      xrange=[6.2,9.8], ytitle='', charsize=1.4, xtitle=''
    im_plothist, mstar[*,ii]-alog10(cat[ii].mu), bin=0.15, /peak, /overplot, $
      /fill, fcolor=im_color('grey80')
    djs_oplot, (ised[ii].mass_50-alog10(cat[ii].mu))*[1,1], !y.crange, $
      line=5, thick=8
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, position=pos2[*,0], yrange=[0,1.05], $
      xrange=[6.2,9.8], ytitle='Probability', charsize=1.4, $
      xtitle=textoidl('log (M_{*}/M'+sunsymbol()+')'), $
      xtickinterval=1, ytickname=replicate(' ',10)

    xr = [-0.04,1] & bin = 0.05
;      xr = [-0.04,0.5] & bin = 0.02
    plot, [0], [0], xsty=5, ysty=5, /noerase, /nodata, yrange=[0,1.05], $
      xrange=xr, position=pos2[*,1], ytitle='', charsize=1.4, $
      xtitle=''
    im_plothist, sfrage[*,ii], bin=bin, /peak, /overplot, $
      /fill, fcolor=im_color('grey80')
    djs_oplot, median(sfrage[*,ii])*[1,1], !y.crange, line=5, thick=8
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=[0,1.05], $
      xrange=xr, position=pos2[*,1], ytitle='', charsize=1.4, $
      xtitle=textoidl('Age_{w} (Gyr)'), ytickname=replicate(' ',10), $
      xtickinterval=0.2

    plot, [0], [0], xsty=5, ysty=5, /noerase, /nodata, yrange=[0,1.1], $
      xrange=[-0.04,1.6], position=pos2[*,2], ytitle='', $
      xtitle='', charsize=1.4
    im_plothist, Z[*,ii]/0.02, bin=0.1, /peak, /overplot, $
      /fill, fcolor=im_color('grey80')
    djs_oplot, ised[ii].Z_50/0.02*[1,1], !y.crange, line=5, thick=8
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=[0,1.1], $
      xrange=[-0.04,1.6], position=pos2[*,2], ytitle='', $
      xtitle='Z/Z'+sunsymbol(), ytickname=replicate(' ',10), charsize=1.4, $
      xtickinterval=0.5

    xr = [-0.03,5.05] & yr = [0,1.1] & xtic = 1
;      xr = [-0.03,2.05] & yr = [0,1.1] & xtic = 0.5
    plot, [0], [0], xsty=5, ysty=5, /noerase, /nodata, yrange=yr, $
      xrange=xr, position=pos2[*,3], ytitle='', $
      xtitle=textoidl('A_{V} (mag)'), $
      ytickname=replicate(' ',10), charsize=1.4, xtickinterval=xtic
    im_plothist, av[*,ii], bin=0.3, /peak, /overplot, $
      /fill, fcolor=im_color('grey80')
    djs_oplot, ised[ii].av_50*[1,1], !y.crange, line=5, thick=8
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=yr, $
      xrange=xr, position=pos2[*,3], ytitle='Probability', $
      xtitle=textoidl('A_{V} (mag)'), $
      ytickname=replicate(' ',10), charsize=1.4, xtickinterval=xtic

    plot, [0], [0], xsty=5, ysty=5, /noerase, /nodata, yrange=[0,1.1], $
      xrange=[-0.5,2.5], position=pos2[*,4], ytitle='', $
      xtitle=textoidl('SFR (M'+sunsymbol()+' yr^{-1})'), $
      ytickname=replicate(' ',10), charsize=1.4, xtickinterval=2
    im_plothist, 10^(sfr0[*,ii]-alog10(cat[ii].mu)), bin=0.2, /peak, /overplot, $
      /fill, fcolor=im_color('grey80')
    djs_oplot, 10^(ised[ii].sfr_50-alog10(cat[ii].mu))*[1,1], !y.crange, line=5, thick=8
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=[0,1.1], $
      xrange=[-0.5,2.5], position=pos2[*,4], ytitle='', $
      xtitle=textoidl('SFR (M'+sunsymbol()+' yr^{-1})'), $
      ytickname=replicate(' ',10), charsize=1.4, xtickinterval=2

    plot, [0], [0], xsty=5, ysty=5, /noerase, /nodata, yrange=[0,1.1], $
      xrange=[0.0,10.0], position=pos2[*,5], ytitle='', $
      ytickname=replicate(' ',10), charsize=1.4
    im_plothist, 10^ssfr[*,ii], bin=0.6, /peak, /overplot, $
      /fill, fcolor=im_color('grey80')
    djs_oplot, median(10^ssfr[*,ii])*[1,1], !y.crange, line=5, thick=8
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=[0,1.1], $
      xrange=[0.0,10.0], position=pos2[*,5], ytitle='', $
      xtitle=textoidl('sSFR (Gyr^{-1})'), $
      ytickname=replicate(' ',10), charsize=1.4, xtickinterval=2

    im_plotconfig, /psclose, psfile=psfile, /pdf, /pskeep

stop


return
end
