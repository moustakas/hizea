function edit_hires_09aug12, rawhires
; jm09aug11ucsd - 

    hires = rawhires
    hires.obj = strtrim(hires.obj,2)
    hires.obj = repstr(hires.obj,' ','_')

    hires.type = strtrim(hires.type,2)

; junk exposures
    crap = where(((hires.frame ge 100) and (hires.frame le 101)) or $
      ((hires.frame ge 115) and (hires.frame le 119)) or $
      ((hires.frame ge 135) and (hires.frame le 146)))
    hires[crap].flg_anly = 0

; BD+33 gets misidentified as an object
    std = where(strmatch(hires.obj,'*BD+33*'),nstd)
    hires[std].type = 'STD'
    
; only keep objects that we need    
    realgood = where(hires.flg_anly eq 1)
    struct_print, struct_trimtags(hires[realgood],sel=['frame',$
      'img_root','obj','colbin','rowbin','type','slit','obj_id',$
      'flg_anly'])

return, hires[realgood]
end
