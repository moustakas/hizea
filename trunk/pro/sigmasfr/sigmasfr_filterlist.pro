function sigmasfr_filterlist
; jm12mar19ucsd 
    filterlist = [galex_filterlist(),sdss_filterlist(),$
      irac_filterlist(/warm),wise_filterlist()]
return, filterlist
end
