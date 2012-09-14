pro hizea_reduce_night, datapath, setup, chip=chip, clobber=clobber, $
  flat=flat, arc=arc, slitflat=slitflat, proc=proc, dotrace=dotrace, $
  skysub=skysub, extract=extract, combine=combine, calibrate=calibrate, $
  dostandards=dostandards, makesens=makesens
; jm09sepucsd - wrapper to reduce one night of hizea/hires data 
    
    if (n_elements(datapath) eq 0L) then datapath = './'
    if (n_elements(setup) eq 0L) then setup = 1 ; default setup
    if (n_elements(chip) eq 0L) then chip = 1   ; [1=blue, 2=green, 3=red]
    
    if (file_test(datapath+'hires.fits',/regular) eq 0L) then begin
       splog, 'HIRES structure not found'
       return
    endif
    hires = hires_ar(datapath+'hires.fits')

    pushd, datapath

; ##################################################
; grab all the objects, and cross-match HIRES against the parent
; catalog of observed objects

; build a hard-coded catalog of all possible objects and cross-match 
;   info = replicate({galaxy: '', ra: 0.0D, dec: 0.0D, g: 0.0D, z: 0.0},5)
;   info.galaxy = ['J1558+3957','J1658+2354','J2116-0634',$
;     'J2140+1209','J2256+1542']
;   info.ra = im_hms2dec(['15:58:00','16:58:00','21:16:25','21:40:00','22:56:00'])*15.0D
;   info.dec = im_hms2dec(['+39:57:00','+23:54:00','-06:34:45','+12:09:15','+15:42:00'])
;
;   info.z = [0.402242,0.498282,0.728,0.751,0.727]
;   info.g = [19.19,19.23,20.03,20.18,20.14]
    
    nobj = 0
    objindx = where((hires.setup eq setup) and (hires.type eq 'OBJ'),nobjindx)
    if (nobjindx ne 0) then begin
       objindx = objindx[uniq(hires[objindx].obj_id,sort(hires[objindx].obj_id))]
;      objindx = objindx[5]
       obj = hires[objindx].obj_id
       nobj = n_elements(obj)
    endif 

; identify the standards and build an INFO structure    
    stdindx = where((hires.setup EQ setup AND hires.type EQ 'STD'),nstd)
    if (nstd gt 0L) then begin
       std = hires[stdindx].obj_id
       stdinfo = struct_trimtags(hires[stdindx],select=['obj_id',$
         'obj','ra','dec','exp','am','img_root'])
;      morestdinfo = hizea_find_standard(stdinfo)
;      stdinfo = struct_addtags(temporary(stdinfo),morestdinfo)
       struct_print, stdinfo
    endif else splog, 'No standard stars observed'

;   niceprint, info.img_root, hires[objindx].img_root
;   niceprint, hires[objindx].frame, hires[objindx].obj, hires[objindx].obj_id

; ##################################################
; process the flats and optionally verify

    if keyword_set(flat) then begin
; (1) generate a stacked, normalized milky flat; (2) combine the trace
; flats for order and slit tracing; (3) generate a smooth model of the
; order curvature using the trace flat previously created
       hires_allflat, hires, setup, chip, /nogain, clobber=clobber
;      spawn, 'gv QA/Flats01/qa_trcflt_01R.ps.gz &'
;      xatv, 'Flats/Flat_R_01_M.fits' 
;      xatv, 'Flats/Flat_R_01_T.fits' 
;      xatv, 'Flats/Flat_R_01_T.fits' 
;      hires_chktrcflat, hires, setup, chip, /nostop, /fit
    endif    

; ##################################################
; process the arcs
    if keyword_set(arc) then begin
       hires_allarc, hires, setup, chip, fits=datapath+'hires.fits', $
         clobber=clobber, chk=chk
    endif

; ##################################################
; build the slit profile
    if keyword_set(slitflat) then begin
       hires_slitflat, hires, setup, chip, clobber=clobber, chk=chk
    endif

; ##################################################
; reduce the standards
    if keyword_set(dostandards) then begin
       if keyword_set(proc) then begin
          for kk = 0L, nstd-1L do hires_proc, hires, setup=setup, $
            obj=std[kk], chip=chip, clobber=clobber, /std
       endif
       if keyword_set(dotrace) then begin
          for kk = 0L, nstd-1L do hires_fntobj, hires, setup, $
            std[kk], chip, chk=chk, clobber=clobber, /std
       endif
       if keyword_set(skysub) then begin
          for kk = 0L, nstd-1L do hires_skysub, hires, setup, $
            stdindx[kk], chip, chk=chk, clobber=clobber, /std ; note STDINDX!
       endif
       if keyword_set(extract) then begin
          for kk = 0L, nstd-1L do hires_box, hires, setup, $
            std[kk], chip, chk=chk, ochk=ochk, reschk=reschk, /std, $
            /skipskysub, nohelio=0, novac=0, clobber=clobber
       endif
    endif

; ##################################################
; build the sensitivity function
    if keyword_set(makesens) then begin
       for kk = 0L, nstd-1L do $
         hires_calibstd, hires, stdindx[kk], $ ; note STDINDX!
         esofil=strtrim(stdinfo[kk].esofil,2) 
    endif 

; ##################################################
; initialize the inverse variance map, overscan-subtract, and divide
; by the flat-field; the output is written in /Final
    if (keyword_set(dostandards) eq 0) and keyword_set(proc) and $
      (nobj ne 0) then begin
       for jj = 0, nobj-1 do begin
          for kk = 0, n_elements(chip)-1 do begin
             hires_proc, hires, setup=setup, obj=obj[jj], $
               chip=chip[kk], clobber=clobber
          endfor
       endfor
    endif
    
; ##################################################
; find and trace the object in the slit; OBJAPER is the aperture to
; mask for sky subtraction (26 unbinned pixels when /STD, and 20
; pixels otherwise); FWIDTH is the fraction of the slit width to use
; when tracing (default 0.25 = 1/4 slit width)
    if (keyword_set(dostandards) eq 0) and keyword_set(dotrace) and $
      (nobj ne 0) then begin
       for jj = 0, nobj-1 do begin
          for kk = 0, n_elements(chip)-1 do begin
;; do some filename juggling to deal with the NOCLOB keyword
;             these = where(obj[jj] eq hires.obj_id)
;             objfil = hires_getfil('obj_fil',/name,frame=hires[these].frame,chip=chip[kk])
             hires_fntobj, hires, setup, obj[jj], chip[kk], $
               objaper=objaper, fwidth=fwidth, chk=chk, noclob=(keyword_set(clobber) eq 0)
          endfor
       endfor  
    endif

; ##################################################
; sky-subtract
    if (keyword_set(dostandards) eq 0) and keyword_set(skysub) and $
      (nobj ne 0) then begin
       for jj = 0, nobj-1 do begin
          for kk = 0, n_elements(chip)-1 do begin
             hires_skysub, hires, setup, obj[jj], chip[kk], chk=chk, clobber=clobber
          endfor
       endfor
    endif

; ##################################################
; extract 1D spectra
    if (keyword_set(dostandards) eq 0) and keyword_set(extract) and $
      (nobj ne 0) then begin
       for jj = 0, nobj-1 do begin
          for kk = 0, n_elements(chip)-1 do begin
;            reschk = 1 & chk = 1
             hires_extract, hires, setup, obj[jj], chip[kk], $
               reschk=reschk, chk=chk, /novac, /modelout, clobber=clobber
          endfor
       endfor
    endif

;; ##################################################
;; flux-calibrate the objects and standard stars
;    if keyword_set(calibrate) then begin
;       if (nstd gt 0L) then begin
;          fluxfil = hires_getfil('sens_fil',setup,chip=chip,$
;            subfil=hires[stdindx[0]].img_root,/name)
;          for kk = 0L, nstd-1L do hires_flux, hires, setup, $
;            stdindx[kk], chip, fluxfil=fluxfil, /std
;          for jj = 0L, nobj-1L do hires_flux, hires, setup, $
;            obj[jj], chip, fluxfil=fluxfil, /boxcar
;       endif else begin
;          splog, 'Unable to flux-calibrate! No standards observed'
;          return
;       endelse
;    endif

; ##################################################
; combine multiple exposures of the same object, which must be run on
; even a single object
   if keyword_set(combine) then begin
;     message, 'Not recommended!!'
      for jj = 0, nobj-1 do begin
         for kk = 0, n_elements(chip)-1 do begin
            hires_combspec, hires, setup, obj[jj], chip[kk], /noflux;, /mchk, /chk, /achk, /noflux
         endfor
      endfor
   endif
   
   popd

return
end
