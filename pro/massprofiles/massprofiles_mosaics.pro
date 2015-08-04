pro massprofiles_mosaics, hst=hst, sdss=sdss
; jm15jul17siena - build the mosaics using scamp + swarp

    massdir = massprofiles_path()
    topdir = massprofiles_path(/hst)
    archdir = topdir+'archive/'
    sdssdir = topdir+'sdss_mosaics/'

    filt = ['F475W','F814W','F160W']
    nfilt = n_elements(filt)

    cat = mrdfits(massdir+'massprofiles_photometry.fits.gz',1)
    ngal = n_elements(cat)

; --------------------------------------------------
; parse the HST data
    if keyword_set(hst) then begin
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
    if keyword_set(sdss) then begin
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
    
;   allgaldir = file_search(archdir+'*',count=ngal)
;   allgal = file_basename(allgaldir)
;   for gg = 0, ngal-1 do begin
;      for ff = 0, nfilt-1 do begin
;         exp = file_search(allgaldir[gg]+'/'+filt[ff]+'/*_flt.fits',$
;           count=nexp)
;         for ee = 0, nexp-1 do file_link, repstr(exp[ee],topdir,''), $
;           topdir+allgal[gg]+'_'+filt[ff]+'_'+strtrim(ee,2)+'.fits'
;      endfor
;   endfor

stop

return
end


