function massprofiles_filterlist, allbands=allbands
; jm15jul06siena
    if keyword_set(allbands) then begin
       filterlist = [galex_filterlist(),sdss_filterlist(),$
         irac_filterlist(/warm),wise_filterlist(shortwave=shortwise)]
    endif else begin
       filterlist = 'clash_wfc3_'+['f390w','f475w','f814w','f160w']+'.par'
    endelse
return, filterlist
end
