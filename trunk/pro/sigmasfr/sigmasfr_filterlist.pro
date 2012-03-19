function sigmasfr_filterlist
; jm12mar19ucsd 
    filterlist = [galex_filterlist(),sdss_filterlist(),wise_filterlist()]
return, filterlist
end
