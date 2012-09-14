function edit_hires_10may20, rawhires
; am10may20ucsd - 

    hires = rawhires
    hires.obj = strtrim(hires.obj,2)
    hires.obj = repstr(hires.obj,' ','_')

    hires.type = strtrim(hires.type,2)

; junk exposures; ignore setups 1.35, 0.9, and 0.8
    crap = where($
      ((hires.frame ge 1) and (hires.frame le 2)) or $   ; test bias
      ((hires.frame ge 3) and (hires.frame le 6)) or $   ; arcs, xdang=1.35
      ((hires.frame ge 11) and (hires.frame le 18)) or $ ; arcs, xdang=0.8, 0.9
      ((hires.frame ge 19) and (hires.frame le 33)) or $ ; flats
      ((hires.frame ge 106) and (hires.frame le 107)) or $ ; feige
      ((hires.frame eq 39))) ; test
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
