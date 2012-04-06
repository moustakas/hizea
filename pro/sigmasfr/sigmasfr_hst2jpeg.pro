pro sigmasfr_hst2jpeg
; jm12apr05ucsd - build jpeg cutouts of the full sample 
    
    path = sigmasfr_path()
    hstpath = sigmasfr_path(/multidrizzle)
    
    hst = rsex(path+'hst_sample.dat')
    lir = mrdfits(path+'sigmasfr_lir.fits.gz',1)
    sigmasfr = alog10(lir.sfr_chary/(!pi*hst.r_e^2)/2.0)
    ngal = n_elements(hst)

    diam_kpc = 30.0 
    
    for ii = 0, ngal-1 do begin
       gal = strlowcase(strtrim(hst[ii].galaxy,2))
       file = file_search(hstpath+'cycle1?/'+gal+'/'+gal+'_drz_sci.fits')
       im = mrdfits(file,0,hdr)

; get the cutout       
       extast, hdr, astr
       sz = size(im,/dim)
       pixscale = sqrt(abs(determ(astr.cd)))*3600D ; [arcsec/pixel]
       diam = round(diam_kpc/dangular(hst[ii].z,/kpc)*206265D/pixscale) ; [pixel]

       ad2xy, hst[ii].ra, hst[ii].dec, astr, xx, yy
       x0 = long(xx-diam/2.0)>0L
       x1 = long(xx+diam/2.0)<(sz[0]-1)
       y0 = long(yy-diam/2.0)>0
       y1 = long(yy+diam/2.0)<(sz[1]-1)
       splog, diam, x0, x1, y0, y1
       
       hextract, im, hdr, newim, newhdr, x0, x1, y0, y1

; make the JPEG
       outfile = path+'jpeg/'+gal+'.png'
       mm = weighted_quantile(newim,quant=[0.0,0.99])
       outim = im_imgscl(newim,/neg,/sqrroot)

       delvarx, pos
       set_plot, 'Z'
       plotimage, outim, /noaxes, /preserve_aspect, position=pos, /norm
       legend, hst[ii].galaxy, /left, /top, box=0, charsize=2.0, $
         charthick=2.0, textcolor=im_color('black',255), $
         margin=0, position=[pos[0]+0.0,pos[3]-0.05], /norm
       legend, textoidl(string(hst[ii].vout,format='(G0)')+' km s^{-1}'), $
         /right, margin=0, /bottom, box=0, charsize=2.0, charthick=2.0, $
         textcolor=im_color('black',255)
;      legend, textoidl([string(10.0^sigmasfr[ii],format='(G0)')+' M'+sunsymbol()+' yr^{-1} kpc^{-2}',$
;        string(hst[ii].vout,format='(G0)')+' km s^{-1}']), $
;        /right, margin=0, /bottom, box=0, charsize=2.2, charthick=2.0, $
;        textcolor=im_color('black',255)
       x0 = fix((convert_coord(pos[0:1],/normal,/to_device))[0])
       nx = fix((convert_coord(pos[2:3],/normal,/to_device))[0])-X0
       y0 = fix((convert_coord(pos[0:1],/normal,/to_device))[1])
       ny = fix((convert_coord(pos[2:3],/normal,/to_device))[1])-y0
       img = tvrd(x0,y0,nx,ny)
       tvlct, rr, gg, bb, /get
       write_png, outfile, img, rr, gg, bb
       set_plot, 'X'

;      outim = asinhscl(newim,/negative,min=mm[0],max=mm[1],beta=0)
;      outim = asinhscl(newim,/negative,min=mm[0],max=mm[1],beta=1.0)
;      write_jpeg, outfile, outim
    endfor

return
end
    
