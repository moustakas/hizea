function edit_hires_12may13, rawhires
; jm12may13ucsd

    hires = rawhires
    hires.obj = strtrim(hires.obj,2)
    hires.obj = repstr(hires.obj,' ','_')

    hires.type = strtrim(hires.type,2)

; junk exposures
    crap = where($
      ((hires.frame ge 28) and (hires.frame le 30)) or $
      ((hires.frame ge 32) and (hires.frame le 36)) or $
      ((hires.frame ge 0118) and (hires.frame le 0119)))
;      (hires.frame eq 0132))
    hires[crap].flg_anly = 0

; ignore pixel flats for now
    toss = where(strtrim(hires.type,2) eq 'PFLT')
    hires[toss].flg_anly = 0
    
; Feige34 gets misidentified as an object
    std = where(strmatch(hires.obj,'*Feige*',/fold) or $
      strmatch(hires.obj,'*BD+28*',/fold),nstd)
    hires[std].type = 'STD'
    
; only keep objects that we need    
    realgood = where(hires.flg_anly eq 1)
    struct_print, struct_trimtags(hires[realgood],sel=['frame',$
      'img_root','obj','colbin','rowbin','type','slit','obj_id',$
      'flg_anly'])

return, hires[realgood]
end
