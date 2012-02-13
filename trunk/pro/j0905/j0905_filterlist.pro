function j0905_filterlist
; jm12feb08ucsd
    filterlist = [galex_filterlist(),sdss_filterlist(),(irac_filterlist())[0:1],$
      'amd_'+['4225','4725','5225','5725','6225','6725','7225','7725']+'.par']
return, filterlist
end
