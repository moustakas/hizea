pro massprofiles_aplus_cutouts, get_cutouts=get_cutouts, trilogy=trilogy
; jm15jul17siena - get cutouts of each galaxy from the APLUS mosaics and
; generate a color image for each using Trilogy

    if n_elements(scale) eq 0 then scale = '65mas'
    
    massdir = massprofiles_path()
    hstdir = massprofiles_path(/hst)

    aplusdir = hstdir+'aplus/'
    cutoutdir = hstdir+'cutouts/'
    colordir = hstdir+'color/'

    filt = ['F475W','F814W','F160W']
    nfilt = n_elements(filt)

    scale = ['25mas','50mas']
    nscale = n_elements(scale)
    
    cat = mrdfits(massdir+'massprofiles_photometry.fits.gz',1)
    ngal = n_elements(cat)

; get cutouts
    if keyword_set(get_cutouts) then begin
       stamp = 350
       for ii = 0, ngal-1 do begin
          gal = strlowcase(strmid(cat[ii].galaxy,0,5)) ; Wei's galaxy name
          for is = 0, nscale-1 do begin
             prefix = strupcase(gal)+'-'+scale[is]
             galdir = aplusdir+prefix
             if file_test(galdir,/dir) eq 0 then begin
                splog, 'Mosaic for '+galdir+' does not (yet) exist.'
             endif else begin
                file_mkdir, hstdir+'cutouts/'+prefix
                for ib = 0, nfilt-1 do begin
                   imfile = galdir+'/'+gal+'_'+strlowcase(filt[ib])+'_sci.fits'
                   rmsfile = repstr(imfile,'_sci.fits','_RMS.fits')
;                  weightfile = repstr(imfile,'_sci.fits','_sci_weight.fits')
                   outfile = hstdir+'cutouts/'+prefix+'/'+prefix+'-'+filt[ib]+'.fits'
                   
                   im = mrdfits(imfile,0,hdr,/silent)
                   rms = mrdfits(rmsfile,0,rmshdr,/silent)
;                  weight = mrdfits(weightfile,0,weighthdr,/silent)
                   
                   extast, hdr, astr
                   ad2xy, cat[ii].ra, cat[ii].dec, astr, xx, yy
                   
                   hextract, im, hdr, newim, newhdr, xx-stamp, xx+stamp, $
                     yy-stamp, yy+stamp, /silent
                   hextract, rms, rmshdr, newrms, newrmshdr, $
                     xx-stamp, xx+stamp, yy-stamp, yy+stamp, /silent
;                  hextract, weight, weighthdr, newweight, newweighthdr, $
;                    xx-stamp, xx+stamp, yy-stamp, yy+stamp, /silent

                   zpt = sxpar(hdr,'ZEROPT')
                   exptime = sxpar(hdr,'EXPTIME')
                   print, filt[ib], zpt, exptime
                   
                   factor = 1D9*10^(-0.4*zpt) ;    [counts-->nanomaggies]
;                  factor = 10^(-0.4*(zpt-31.4)) ;    [counts-->nJy]
                   splog, 'Writing '+outfile
                   mwrfits, float(newim*factor), outfile, newhdr, /create
                   
                   ivar = 1D/(newrms^2 + (newrms eq 0))*(newrms ne 0)
                   mwrfits, float(ivar*1/factor^2), $
                     repstr(outfile,'.fits','-ivar.fits'), $
                     newrmshdr, /create
;                  mwrfits, float(newweight*1/factor^2), $
;                    repstr(outfile,'.fits','-ivar.fits'), $
;                    newweighthdr, /create
                   
;   ;    also get the PSF
;                  psffile = aplusdir+prefix+'-'+scale+'/'+strlowcase(prefix)+'_'+$
;                    strlowcase(filt[ib])+'.psf'
;                  psf = mrdfits(psffile,1,psfhdr)
;                  
;                  outpsffile = aplusdir+'cutouts/'+gal[ii]+'-'+scale+'-'+filt[ib]+'.psf'
;                  mwrfits, psf, outpsffile, psfhdr, /create
                endfor
             endelse
          endfor 
       endfor
    endif

; call trilogy
    if keyword_set(trilogy) then begin
       for ii = 0, ngal-1 do begin
          for is = 0, nscale-1 do begin
             prefix = strupcase(strmid(cat[ii].galaxy,0,5))+'-'+scale[is]
             cutdir = hstdir+'cutouts/'+prefix
             if file_test(cutdir,/dir) ne 0 then begin
                openw, lun, colordir+'trilogy.in', /get_lun
                printf, lun, 'B'
                printf, lun, prefix+'-'+filt[0]
                printf, lun, ''
                printf, lun, 'G'
                printf, lun, prefix+'-'+filt[1]
                printf, lun, ''
                printf, lun, 'R'
                printf, lun, prefix+'-'+filt[2]
                printf, lun, ''
                printf, lun, 'indir '+cutdir
                printf, lun, 'outdir '+colordir
                printf, lun, 'outname '+prefix
                printf, lun, 'noiselum 0.3'
                printf, lun, 'satpercent  0.001'
                printf, lun, 'legend 0'
                printf, lun, 'deletetests 1'
;               printf, lun, 'samplesize 100'
                printf, lun, 'show 0'
                printf, lun, 'testfirst 0'
                free_lun, lun
                spawn, 'trilogy '+colordir+'trilogy.in'
                file_delete, colordir+'trilogy.in'
                file_delete, colordir+'levels.txt'
                file_delete, colordir+prefix+'_filters.txt'
                file_delete, 'trilogyfilterlog.txt'
             endif
          endfor 
       endfor 
    endif 
    
return
end


