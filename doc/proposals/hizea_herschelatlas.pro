pro hizea_herschelatlas

; 15 degrees by 10 degrees centred on RA=199.5, Dec=29 and rotated by 8
; 26.7 degrees in RA by 6 degrees in Dec, centred on RA=351.3, Dec=-32.8
; 17 degrees by 6 degrees centred on RA=36.7, Dec=-30.7

    jj = mrdfits(getenv('HIZEA_DATA')+'/sdss/hizea_sdss_photo_dr7.fits.gz',1)
    hst = rsex(getenv('HIZEA_DATA')+'/sigmasfr/hst_sample.dat')
    
    djs_plot, [0], [0], xrange=[0,360]/15D, yrange=[-20,90], xsty=1, ysty=1
    djs_oplot, jj.ra/15D, jj.dec, psym=symcat(6,thick=1), symsize=0.5
    djs_oplot, hst.ra/15D, hst.dec, psym=symcat(6), symsize=0.5, color='red'

    im_oplot_box, 15.0/15D, 10.0, 8.0, xoffset=(199.5-15/2.0)/15D, $
      yoffset=29.0-10.0/2.0, thick=4, color='orange'
;   im_oplot_box, 26.7, 6.0, 0.0, xoffset=351.3, yoffset=-32.8, thick=6, color='orange'
;   im_oplot_box, 17.0, 6.0, 0.0, xoffset=36.7, yoffset=-30.7, thick=6, color='orange'

    
return
end
