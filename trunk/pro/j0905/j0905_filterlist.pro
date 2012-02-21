function j0905_filterlist, nomediumband=nomediumband
; jm12feb08ucsd
    filterlist = [galex_filterlist(),sdss_filterlist(),(irac_filterlist())[0:1],$
      wise_filterlist(),$
      'amd_'+['4225','4725','5225','5725','6225','6725','7225','7725']+'_v2.par']
    if keyword_set(nomediumband) then filterlist = $
      filterlist[where(strmatch(filterlist,'*amd*',/fold) eq 0)]
return, filterlist
end
