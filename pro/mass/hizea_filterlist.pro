function hizea_filterlist
; jm10dec20ucsd
    filterlist = [galex_filterlist(),sdss_filterlist(),irac_filterlist()]
return, filterlist
end
