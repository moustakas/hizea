pro massprofiles_sort_archive, clobber=clobber
; jm15jul17siena - sort/organize the data downloaded from the archive plus the
; SDSS mosaics by galaxy; this is the key preprocessing step before Wei runs the
; APLUS pipeline on the data

    if n_elements(scale) eq 0 then scale = '65mas'
    
    massdir = massprofiles_path()
    hstdir = massprofiles_path(/hst)

    aplusdir = hstdir+'aplus/'
    archdir = hstdir+'archive/'
    sorteddir = hstdir+'archive_sorted/'
    sdssdir = hstdir+'sdss_mosaics/'

    filt = ['F475W','F814W','F160W']
    nfilt = n_elements(filt)

    cat = mrdfits(massdir+'massprofiles_photometry.fits.gz',1)
    ngal = n_elements(cat)

    gal = ['J1341-0321']
    
    allfits = file_search(archdir+'*_flt.fits',count=nall)
    info = replicate({file: '', rootname: '', galaxy: '', $
      ra: 0D, dec: 0D, filter: ''},nall)
    
; forage the FITS headers and then group by ra,dec,filter
    for ii = 0, nall-1 do begin
       hdr = headfits(allfits[ii])
       info[ii].file = file_basename(allfits[ii])
       info[ii].rootname = strtrim(sxpar(hdr,'ROOTNAME'),2)
       info[ii].galaxy = strtrim(sxpar(hdr,'TARGNAME'),2)
       info[ii].ra = sxpar(hdr,'RA_TARG')
       info[ii].dec = sxpar(hdr,'DEC_TARG')
       info[ii].filter = strtrim(sxpar(hdr,'FILTER'),2)
    endfor
    info = info[sort(info.ra)]
    outfile = hstdir+'archive-info.fits'
    if im_file_test(outfile,clobber=clobber) eq 0 then begin
       splog, 'Writing '+outfile
       mwrfits, info, outfile, /create
    endif
    
    for ic = 0, ngal-1 do begin
       spherematch, info.ra, info.dec, cat[ic].ra, cat[ic].dec, $
         5D/3600.0, m1, m2, maxmatch=0
       gal = strtrim(cat[ic].galaxy,2)
       file_mkdir, sorteddir+gal
       for ff = 0, nfilt-1 do begin
          these = where(filt[ff] eq strtrim(info[m1].filter,2),nexp)
          for nn = 0, nexp-1 do begin
             outfile = sorteddir+gal+'/'+file_basename(strtrim(info[m1[these[nn]]].file,2))
;            outfile = hstdir+gal+'/'+gal+'-'+filt[ff]+'-'+strtrim(nn,2)+'.fits'
;            splog, strtrim(info[m1[these[nn]]].file,2)+' --> '+outfile
             if im_file_test(outfile,clobber=clobber) eq 0 then begin
                splog, 'Writing '+outfile
                file_copy, archdir+strtrim(info[m1[these[nn]]].file,2), outfile, /over
             endif
          endfor
       endfor 
       print
    endfor
    
; parse the SDSS mosaics
    allfits = file_search(sdssdir+'*-r.fits',count=nall)
    info = replicate({file: '', galaxy: '', ra: 0D, $
      dec: 0D, filter: 'r'},nall)

; forage the FITS headers and then group by ra,dec,filter
    for ii = 0, nall-1 do begin
       hdr = headfits(allfits[ii])
       extast, hdr, astr
       info[ii].file = allfits[ii]
;      info[ii].galaxy = sxpar(hdr,'TARGNAME')
       info[ii].ra = astr.crval[0]
       info[ii].dec = astr.crval[1]
    endfor

    spherematch, info.ra, info.dec, cat.ra, cat.dec, $
      5D/3600.0, m1, m2
    if n_elements(m1) ne n_elements(m2) then message, 'Problem here!'
    
    for ii = 0, ngal-1 do begin
       gal = strtrim(cat[m2[ii]].galaxy,2)
       outfile = sorteddir+gal+'/sdss-'+gal+'.fits'
       splog, 'Writing '+outfile
       if im_file_test(outfile,clobber=clobber) eq 0 then begin
          splog, 'Writing '+outfile
          file_copy, strtrim(info[m1[ii]].file,2), outfile, /over
       endif
    endfor
    
return
end


