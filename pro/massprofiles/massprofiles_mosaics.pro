pro massprofiles_mosaics
; jm15jul17siena - build the mosaics using scamp + swarp

    topdir = massprofiles_path(/hst)
    rawdir = topdir+'rawdata/'

    filt = ['F475W','F814W','F160W']
    nfilt = n_elements(filt)

    allgaldir = file_search(rawdir+'*',count=ngal)
    allgal = file_basename(allgaldir)
    for gg = 0, ngal-1 do begin
       for ff = 0, nfilt-1 do begin
          exp = file_search(allgaldir[gg]+'/'+filt[ff]+'/*_flt.fits',$
            count=nexp)
          for ee = 0, nexp-1 do file_link, repstr(exp[ee],topdir,''), $
            topdir+allgal[gg]+'_'+filt[ff]+'_'+strtrim(ee,2)+'.fits'
       endfor
    endfor

stop

return
end


