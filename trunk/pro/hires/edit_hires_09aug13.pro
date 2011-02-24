function edit_hires_09aug13, rawhires
; jm09sep01ucsd - 

    hires = rawhires
    hires.obj = strtrim(hires.obj,2)
    hires.obj = repstr(hires.obj,' ','_')

    hires.type = strtrim(hires.type,2)

; don't analyze the pixel flats here
    pflat = where(hires.type eq 'PFLT')
    hires[pflat].flg_anly = 0
    
; remove the Setup A arcs and quartz flats since no science objects
; were observed using this configuration
    setupa = where(((hires.frame ge 1102) and (hires.frame le 1106)) or $
      ((hires.frame ge 1117) and (hires.frame le 1121)))
    hires[setupa].flg_anly = 0
    
; only keep objects that we need    
    realgood = where(hires.flg_anly eq 1)
    struct_print, struct_trimtags(hires[realgood],sel=['frame',$
      'img_root','obj','colbin','rowbin','type','slit','obj_id',$
      'flg_anly'])

return, hires[realgood]
end
