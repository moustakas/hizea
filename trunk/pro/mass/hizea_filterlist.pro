function hizea_filterlist
; jm10dec20ucsd
    filterlist = [galex_filterlist(),sdss_filterlist(),(irac_filterlist())[0:1]]
return, filterlist
end
