function hizea_filterlist, shortwise=shortwise
; jm10dec20ucsd
; jm15nov16siena - updated to include WISE
    filterlist = [galex_filterlist(),sdss_filterlist(),$
      irac_filterlist(/warm),wise_filterlist(shortwave=shortwise)]
return, filterlist
end
