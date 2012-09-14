function edit_hires_12may14, rawhires
; jm12may14ucsd

    hires = rawhires
    hires.obj = strtrim(hires.obj,2)
    hires.obj = repstr(hires.obj,' ','_')

    hires.type = strtrim(hires.type,2)

; junk exposures
    crap = where($
      ((hires.frame ge 200) and (hires.frame le 202)) or $ ; arc covers closed
      (hires.frame eq 200) or $
      (hires.frame eq 292) or $ ; feige34 at xdangle=1.1
      ((hires.frame ge 212) and (hires.frame le 220)) or $ ; calibrations at xdangle=1.1
      ((hires.frame ge 261) and (hires.frame le 290)))
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
