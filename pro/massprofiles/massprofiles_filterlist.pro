function massprofiles_filterlist, shortwise=shortwise
; jm15jul06siena
    filterlist = [galex_filterlist(),sdss_filterlist(),$
      irac_filterlist(/warm),wise_filterlist(shortwave=shortwise)]
return, filterlist
end
