pro massprofiles_apphot
; jm15aug12siena - do aperture photometry on the APLUS cutouts and build the
; full SED

    massdir = massprofiles_path()
    topdir = massprofiles_path(/hst)
    aplusdir = topdir+'aplus/cutouts/'

    cat = mrdfits(massdir+'massprofiles_photometry.fits.gz',1)
    ngal = n_elements(cat)

    filt = ['F475W','F814W','F160W']
    filterlist = 'clash_wfc3_'+['f475w','f814w','f160w']+'.par'
    weff = k_lambda_eff(filterlist=filterlist)
    nfilt = n_elements(filt)

    allfilt = sigmasfr_filterlist()
    allweff = k_lambda_eff(filterlist=allfilt)

    nrad = 10
    rad = range(10,100,nrad) ; [pixels]
    rad_arcsec = rad*0.030   ; [arcsec]

    gal = ['J1341-0321']
    for ii = 0, ngal-1 do begin
       this = where(strmatch(cat.galaxy,gal[ii]+'*'))

       phot = {rad: rad, rad_arcsec: rad_arcsec, $
         maggies: fltarr(nfilt,nrad), ivarmaggies: fltarr(nfilt,nrad)}
       
       for ff = 0, nfilt-1 do begin
          im = mrdfits(aplusdir+gal[ii]+'-'+filt[ff]+'.fits',0,hdr)
          ivar = mrdfits(aplusdir+gal[ii]+'-'+filt[ff]+'-ivar.fits',0,hdr)
          sz = size(im,/dim)
          xcen = sz[0]/2.0
          ycen = sz[1]/2.0

;         plotimage, im_imgscl(im), /preserve
;         tvcircle, 30, xcen, ycen, /data, color=cgcolor('green')
;         tvcircle, 70, xcen, ycen, /data, color=cgcolor('green')

          phot.maggies[ff,*] = djs_phot(xcen,ycen,rad,$
            skyrad,im,ivar,calg='none',flerr=ferr1)
          phot.ivarmaggies[ff,*] = 1/ferr1^2
       endfor

       psfile = massdir+'qa-'+gal[ii]+'-sed.eps'
       im_plotconfig, 0, pos, psfile=psfile, height=5.0
       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         /xlog, xrange=[0.1,30], yrange=[22.5,14], $
         xtitle='Wavelength (\mu'+'m)', ytitle='AB Magnitude'
       good = where(cat[this].maggies gt 0)
       mag = maggies2mag(cat[this].maggies[good],$
         ivarmaggies=cat[this].ivarmaggies[good],$
         magerr=magerr)
       oploterror, allweff[good]/1D4, mag, magerr, psym=cgsymcat(15), $
         symsize=2, errthick=4
       mag = maggies2mag(phot.maggies[*,6],magerr=magerr)
;        ivarmaggies=phot.ivarmaggies[*,6])
       djs_oplot, weff/1D4, mag, psym=cgsymcat(16), $
         symsize=2, color=cgcolor('orange')
;      oploterror, weff/1D4, mag, magerr, psym=cgsymcat(16), $
;        symsize=2, color=cgcolor('orange'), errthick=4
       im_legend, ['HST/WFC3','GALEX+SDSS+WISE'], /left, /top, $
         box=0, color=['orange','black'], psym=[16,15]
       im_plotconfig, psfile=psfile, /psclose, /png

stop       
       
       psfile = massdir+'qa-'+gal[ii]+'-rad.eps'
       im_plotconfig, 0, pos, psfile=psfile, height=5.0
       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         /xlog, xrange=[0.1,30], yrange=[22.5,14], $
         xtitle='Wavelength (\mu'+'m)', ytitle='AB Magnitude'
       good = where(cat[this].maggies gt 0)
       mag = maggies2mag(cat[this].maggies[good],$
         ivarmaggies=cat[this].ivarmaggies[good],$
         magerr=magerr)
       oploterror, allweff[good]/1D4, mag, magerr, psym=cgsymcat(15), $
         symsize=2, errthick=4
       mag = maggies2mag(phot.maggies[*,6],magerr=magerr,$
         ivarmaggies=phot.ivarmaggies[*,6])
       oploterror, weff/1D4, mag, magerr, psym=cgsymcat(16), $
         symsize=2, color=cgcolor('orange'), errthick=4
       im_plotconfig, psfile=psfile, /psclose, /png
       
       
stop       
       
    endfor
    
return
end
    
