function edit_hires_10may20, rawhires
; am10may20ucsd - 

    hires = rawhires
    hires.obj = strtrim(hires.obj,2)
    hires.obj = repstr(hires.obj,' ','_')
    hires.type = strtrim(hires.type,2)

; Mendez - Removing Pixel Flats and the sort,  I am not sure how they
;   are handled, so for the "quick" reduction, ignoring.
  crap = where((hires.frame ge 39) and (hires.frame le 104), ncrap)
  if ncrap gt 0 then hires[crap].flg_anly=0
  
; Mendez - Removing Bias Frames,  so that no hires.setup==99
  bias = where((hires.frame ge 1) and (hires.frame le 2), nbias)
  if nbias gt 0 then hires[bias].flg_anly=0

; Mendez - Assigning STD Star: 'Feige 34':
;  I did not check to see that this was aassigned, just asigned it.
  std = where(strmatch(hires.obj,'*feige*',/fold),nstd)
  if nstd gt 0 then hires[std].type = 'STD'
  ; if nstd gt 0 then hires[std].flg_anly=0 ; Trying to see if that works

; Mendez - First 30min exp, has an order that wont reduce,  rm for now
;   There one can ignore the bad orders with fix in ~ajmendez/idl/hizea_redux
;   This is only to process a single frame, so that when I rsync 
;   a new file, I can process it quite quickly ~3.6min
  badframe = where(hires.frame eq 108 or $
                   hires.frame eq 109 or $
                   hires.frame eq 110 or $
                   hires.frame eq 111 or $
                   hires.frame eq 112 or $
                   hires.frame eq 113 or $
                   hires.frame eq 114 or $
                   hires.frame eq 115, nbad)
  if nbad gt 0 then hires[badframe].flg_anly=0
  
; only keep objects that we need    
    realgood = where(hires.flg_anly eq 1)
    struct_print, struct_trimtags(hires[realgood],sel=['frame',$
      'img_root','obj','colbin','rowbin','type','slit','obj_id',$
      'flg_anly'])

;Mendez - Turn on stop, to check, then off after any corrections.
  ; stop
return, hires[realgood]
end
