function edit_hires_09dec10, rawhires
; jm09dec10ucsd - 

    hires = rawhires
    hires.obj = strtrim(hires.obj,2)
    hires.obj = repstr(hires.obj,' ','_')

    hires.type = strtrim(hires.type,2)

; remove test bias frames
    crap = where((hires.frame eq 101) or (hires.frame eq 102)) 
    hires[crap].flg_anly = 0

; Feige 34 is misidentified as an OBJECT
    std = where((strmatch(hires.obj,'*feige*',/fold)))
    hires[std].type = 'STD'

; only keep objects that we need    
    realgood = where(hires.flg_anly eq 1)
    struct_print, struct_trimtags(hires[realgood],sel=['frame',$
      'img_root','obj','colbin','rowbin','type','slit','obj_id',$
      'flg_anly'])

return, hires[realgood]
end
