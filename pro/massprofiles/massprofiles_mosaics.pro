pro massprofiles_mosaics, parse_hst=parse_hst, parse_sdss=parse_sdss, $
  get_cutouts=get_cutouts, trilogy=trilogy, build_psf=build_psf, $
  sextractor=sextractor, scale=scale, dual=dual
; jm15jul17siena - process the raw data so that Wei can run APLUS

    if n_elements(scale) eq 0 then scale = '65mas'
    
    massdir = massprofiles_path()
    topdir = massprofiles_path(/hst)
    hizeadir = massprofiles_path(/code)+'etc/'

    aplusdir = topdir+'aplus/'
    archdir = topdir+'archive/'
    sdssdir = topdir+'sdss_mosaics/'

    filt = ['F475W','F814W','F160W']
    nfilt = n_elements(filt)

    cat = mrdfits(massdir+'massprofiles_photometry.fits.gz',1)
    ngal = n_elements(cat)

    gal = ['J1341-0321']
    
; --------------------------------------------------
; parse the HST data
    if keyword_set(parse_hst) then begin
       allfits = file_search(archdir+'*_flt.fits',count=nall)
       info = replicate({file: '', rootname: '', galaxy: '', $
         ra: 0D, dec: 0D, filter: ''},nall)

; forage the FITS headers and then group by ra,dec,filter
       for ii = 0, nall-1 do begin
          hdr = headfits(allfits[ii])
          info[ii].file = allfits[ii]
          info[ii].rootname = sxpar(hdr,'ROOTNAME')
          info[ii].galaxy = sxpar(hdr,'TARGNAME')
          info[ii].ra = sxpar(hdr,'RA_TARG')
          info[ii].dec = sxpar(hdr,'DEC_TARG')
          info[ii].filter = sxpar(hdr,'FILTER')
       endfor
       info = info[sort(info.ra)]
       mwrfits, info, topdir+'archive-info.fits', /create
       
       for ic = 0, ngal-1 do begin
          spherematch, info.ra, info.dec, cat[ic].ra, cat[ic].dec, $
            5D/3600.0, m1, m2, maxmatch=0
          gal = strtrim(cat[ic].galaxy,2)
          file_mkdir, topdir+gal
          for ff = 0, nfilt-1 do begin
             these = where(filt[ff] eq strtrim(info[m1].filter,2),nexp)
             for nn = 0, nexp-1 do begin
                outfile = topdir+gal+'/'+file_basename(strtrim(info[m1[these[nn]]].file,2))
;               outfile = topdir+gal+'/'+gal+'-'+filt[ff]+'-'+strtrim(nn,2)+'.fits'
;               splog, strtrim(info[m1[these[nn]]].file,2)+' --> '+outfile
                splog, 'Writing '+outfile
                file_copy, strtrim(info[m1[these[nn]]].file,2), outfile, /over
             endfor
          endfor 
          print
       endfor
    endif
    
; --------------------------------------------------
; parse the SDSS mosaics
    if keyword_set(parse_sdss) then begin
       allfits = file_search(sdssdir+'*-r.fits',count=nall)
       info = replicate({file: '', galaxy: '', ra: 0D, $
         dec: 0D, filter: 'r'},nall)

; forage the FITS headers and then group by ra,dec,filter
       for ii = 0, nall-1 do begin
          hdr = headfits(allfits[ii])
          extast, hdr, astr
          info[ii].file = allfits[ii]
;         info[ii].galaxy = sxpar(hdr,'TARGNAME')
          info[ii].ra = astr.crval[0]
          info[ii].dec = astr.crval[1]
       endfor

       spherematch, info.ra, info.dec, cat.ra, cat.dec, $
         5D/3600.0, m1, m2
       if n_elements(m1) ne n_elements(m2) then message, 'Problem here!'
       
       for ii = 0, ngal-1 do begin
          gal = strtrim(cat[m2[ii]].galaxy,2)
          outfile = topdir+gal+'/sdss-'+gal+'.fits'
          splog, 'Writing '+outfile
          file_copy, strtrim(info[m1[ii]].file,2), outfile, /over
       endfor
    endif 

; --------------------------------------------------
; run SExtractor and then build the PSF using PSFEx
    if keyword_set(sextractor) then begin
       for ii = 0, n_elements(gal)-1 do begin
          prefix = strlowcase(strmid(gal[ii],0,5))
; next run SE in dual-image mode (or not!)
          detect_imfile = aplusdir+prefix+'-'+scale+'/'+prefix+$
            '_f160w_sci.fits'
          detect_rmsfile = repstr(detect_imfile,'_sci.fits','_RMS.fits')
          for ib = 0, nfilt-1 do begin
             imfile = aplusdir+prefix+'-'+scale+'/'+prefix+'_'+$
               strlowcase(filt[ib])+'_sci.fits'
             rmsfile = repstr(imfile,'_sci.fits','_RMS.fits')
             flagfile = repstr(imfile,'_sci.fits','_FLAG.fits')
             catfile = repstr(imfile,'_sci.fits','.cat')

             if keyword_set(dual) then begin
                cmd = 'sex '+detect_imfile+','+imfile+' '+$
                  '-c '+hizeadir+'massprofiles.se '+$
                  '-PARAMETERS_NAME '+hizeadir+'massprofiles.param '+$
                  '-FILTER_NAME '+hizeadir+'gauss_1.5_3x3.conv '+$
                  '-STARNNW_NAME '+hizeadir+'default.nnw '+$
                  '-WEIGHT_TYPE MAP_RMS,MAP_RMS '+$
                  '-WEIGHT_IMAGE '+detect_rmsfile+','+rmsfile+' '+$
                  '-FLAG_IMAGE '+flagfile+' '+$
                  '-CATALOG_NAME '+catfile
                spawn, cmd, /sh
             endif else begin
                cmd = 'sex '+imfile+' '+$
                  '-c '+hizeadir+'massprofiles.se '+$
                  '-PARAMETERS_NAME '+hizeadir+'massprofiles.param '+$
                  '-FILTER_NAME '+hizeadir+'gauss_1.5_3x3.conv '+$
                  '-STARNNW_NAME '+hizeadir+'default.nnw '+$
                  '-WEIGHT_TYPE MAP_RMS '+$
                  '-WEIGHT_IMAGE '+rmsfile+' '+$
                  '-FLAG_IMAGE '+flagfile+' '+$
                  '-CATALOG_NAME '+catfile
                spawn, cmd, /sh
             endelse
; write a DS9 region file
             cat = mrdfits(catfile,2)
             write_ds9_regionfile, cat.alpha_j2000, cat.delta_j2000, $
               filename='j1341'+strlowcase(filt[ib])+'.reg'
          endfor
       endfor
    endif

    if keyword_set(build_psf) then begin
       for ii = 0, n_elements(gal)-1 do begin
          prefix = strlowcase(strmid(gal[ii],0,5))
          for ib = 0, nfilt-1 do begin
             catfile = aplusdir+prefix+'-'+scale+'/'+prefix+'_'+$
               strlowcase(filt[ib])+'.cat'
             outcat = repstr(catfile,'.cat','-psf.cat')
             sampfile = repstr(catfile,'.cat','-psf-samp.fits')
             residfile = repstr(catfile,'.cat','-psf-resid.fits')
             cmd = 'psfex '+$
               '-c '+hizeadir+'massprofiles.psfex '+$
               '-OUTCAT_NAME '+outcat+' '+$
               '-CHECKIMAGE_NAME samp,resid '+$
;              '-CHECKIMAGE_NAME '+sampfile+','+residfile+' '+$
               catfile
             spawn, cmd, /sh
          endfor
       endfor               
    endif

; --------------------------------------------------
; get cutouts
    if keyword_set(get_cutouts) then begin
       stamp = 300

       for ii = 0, n_elements(gal)-1 do begin
          prefix = strlowcase(strmid(gal[ii],0,5))
          this = where(strmatch(cat.galaxy,gal[ii]+'*'))
          for ib = 0, nfilt-1 do begin
             imfile = aplusdir+prefix+'-'+scale+'/'+strlowcase(prefix)+'_'+$
               strlowcase(filt[ib])+'_sci.fits'
             rmsfile = repstr(imfile,'_sci.fits','_RMS.fits')
;            weightfile = repstr(imfile,'_sci.fits','_sci_weight.fits')
             outfile = aplusdir+'cutouts/'+gal[ii]+'-'+scale+'-'+filt[ib]+'.fits'

             im = mrdfits(imfile,0,hdr,/silent)
             rms = mrdfits(rmsfile,0,rmshdr,/silent)
;            weight = mrdfits(weightfile,0,weighthdr,/silent)

             extast, hdr, astr
             ad2xy, cat[this].ra, cat[this].dec, astr, xx, yy

             hextract, im, hdr, newim, newhdr, xx-stamp, xx+stamp, $
               yy-stamp, yy+stamp, /silent
             hextract, rms, rmshdr, newrms, newrmshdr, $
               xx-stamp, xx+stamp, yy-stamp, yy+stamp, /silent
;            hextract, weight, weighthdr, newweight, newweighthdr, $
;              xx-stamp, xx+stamp, yy-stamp, yy+stamp, /silent

             zpt = sxpar(hdr,'ZEROPT')
             exptime = sxpar(hdr,'EXPTIME')
             print, filt[ib], zpt, exptime

             factor = 1D9*10^(-0.4*zpt) ; [counts-->nanomaggies]
;            factor = 10^(-0.4*(zpt-31.4)) ; [counts-->nJy]
             mwrfits, float(newim*factor), outfile, newhdr, /create

             ivar = 1D/(newrms^2 + (newrms eq 0))*(newrms ne 0)
             mwrfits, float(ivar*1/factor^2), $
               repstr(outfile,'.fits','-ivar.fits'), $
               newrmshdr, /create
;            mwrfits, float(newweight*1/factor^2), $
;              repstr(outfile,'.fits','-ivar.fits'), $
;              newweighthdr, /create

; also get the PSF
             psffile = aplusdir+prefix+'-'+scale+'/'+strlowcase(prefix)+'_'+$
               strlowcase(filt[ib])+'.psf'
             psf = mrdfits(psffile,1,psfhdr)
             
             outpsffile = aplusdir+'cutouts/'+gal[ii]+'-'+scale+'-'+filt[ib]+'.psf'
             mwrfits, psf, outpsffile, psfhdr, /create
          endfor 
       endfor
    endif

; --------------------------------------------------
; call trilogy
    if keyword_set(trilogy) then begin
       for ii = 0, n_elements(gal)-1 do begin
          openw, lun, aplusdir+'cutouts/trilogy.in', /get_lun
          printf, lun, 'B'
          printf, lun, gal[ii]+'-'+scale+'-'+filt[0]
          printf, lun, ''
          printf, lun, 'G'
          printf, lun, gal[ii]+'-'+scale+'-'+filt[1]
          printf, lun, ''
          printf, lun, 'R'
          printf, lun, gal[ii]+'-'+scale+'-'+filt[2]
          printf, lun, ''
          printf, lun, 'indir '+aplusdir+'cutouts'
          printf, lun, 'outdir '+aplusdir+'cutouts'
          printf, lun, 'outname '+gal[ii]+'-'+scale
          printf, lun, 'noiselum 0.3'
          printf, lun, 'satpercent  0.001'
          printf, lun, 'legend 0'
          printf, lun, 'deletetests 1'
;         printf, lun, 'samplesize 100'
          printf, lun, 'show 0'
          printf, lun, 'testfirst 0'
          free_lun, lun
          spawn, 'trilogy '+aplusdir+'cutouts/trilogy.in'
          file_delete, aplusdir+'cutouts/trilogy.in'
          file_delete, aplusdir+'cutouts/levels.txt'
          file_delete, aplusdir+'cutouts/'+gal[ii]+'-'+scale+'_filters.txt'
          file_delete, 'trilogyfilterlog.txt'
       endfor
    endif 

    
return
end


