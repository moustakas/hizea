pro hizea_chandra_plots
; jm13sep04siena - build some plots

    prefix = 'hizea_chandra'
    isedfit_dir = getenv('HIZEA_DATA')+'/chandra/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    params = read_isedfit_paramfile(isedfit_paramfile)

    cat = mrdfits(isedfit_dir+'hizea_chandra_photometry.fits.gz',1)
    isedfit_results = read_isedfit(isedfit_paramfile,/flam,$
      isedfit_dir=isedfit_dir,thissfhgrid=1,/getmodels)
    ngal = n_elements(isedfit_results)
    
; ---------------------------------------------------------------------------    
; photometric best-fits vs spectroscopy
    filterlist = strtrim(params.filterlist,2)
    nfilt = n_elements(filterlist)
    weff = k_lambda_eff(filterlist=filterlist)
    hwhm = weff*0.0

    yscale = 1D17
    wscale = 1D4
    light = im_light(/ang)
    
    yrange = [0,20]
    xrange2 = [1500.0,7000]/wscale
    
    psfile = isedfit_dir+'qa_photspec.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0
    for igal = 0, ngal-1L do begin
       z = isedfit_results[igal].z
       wave = isedfit_results[igal].wave    ; [A]
       flux = isedfit_results[igal].flux    ; [AB]

       xrange1 = xrange2*(1.0+z)
       
       these = where(isedfit_results[igal].ivarmaggies gt 0)
       factor = 10^(-0.4*48.6)*light/weff[these]^2

       best = isedfit_results[igal].bestmaggies[these]*factor
       obs = isedfit_results[igal].maggies[these]*factor
       obserr = 1.0/sqrt(isedfit_results[igal].ivarmaggies[these])*factor

       djs_plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
         xsty=1, ysty=1, xtitle='Wavelength \lambda (\mu'+'m)', $
         ytitle='F_{\lambda} (10^{-17} '+flam_units()+')', position=pos
       im_legend, cat[igal].galaxy, /left, /bottom, box=0

; get the MMT spectrum
       mmt = mrdfits(isedfit_dir+'mmt/'+strtrim(cat[igal].galaxy,2)+'_spec.fits',1,/silent)
       djs_oplot, mmt.wave/wscale, mmt.flux, color='orange'
       
       djs_oplot, wave/wscale, flux*yscale, line=0, color='grey'
       djs_oplot, weff[these]/wscale, best*yscale, $
         psym=symcat(6,thick=6), symsize=2.5
       oploterror, weff[these]/wscale, obs*yscale, obserr*yscale, psym=symcat(16), $
         symsize=1.7, color=im_color('dodger blue'), $
         errcolor=im_color('dodger blue'), errthick=!p.thick

    endfor 
    im_plotconfig, psfile=psfile, /psclose, /pdf
    

stop    
    
; ---------------------------------------------------------------------------    
; SFH plots    
    loadct, 13, /silent
    color = fix((findgen(ngal)*(255-0.0)+0.0)/float(ngal))

    psfile = isedfit_dir+'qa_sfh.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=4.7, $
      xmargin=[1.8,0.4], width=6.3

    xrange = [-5.0,1.0]
;   xrange = [0.1,9.5]
    djs_plot, [0], [0], /nodata, xrange=xrange, yrange=[1E-2,2E3], $
      /ylog, xsty=1, ysty=1, position=pos, xtitle='Time from Burst (Gyr)', $
      ytitle='SFR (M_{'+sunsymbol()+'} yr^{-1})'
    for ii = 0, ngal-1 do begin
       delvarx, aa
       sfh = isedfit_sfh(ised[ii],maxage=9.5,outage=aa,/linear)
       toff = ised[ii].tburst
       before = where(aa lt ised[ii].age,comp=after)
       cgoplot, aa[before]-toff, sfh[before], line=0, thick=6, color=color[ii]
       cgoplot, aa[after]-toff, sfh[after], line=2, thick=2, color=color[ii]
    endfor

    im_plotconfig, psfile=psfile, /psclose, /pdf
    
stop

return
end
    
