function sigmasfr_filterlist, shortwise=shortwise
; jm12mar19ucsd 
    filterlist = [galex_filterlist(),sdss_filterlist(),$
      irac_filterlist(/warm),wise_filterlist(shortwave=shortwise)]
return, filterlist
end
