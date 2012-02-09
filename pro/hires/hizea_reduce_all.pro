pro hizea_reduce_all, night, setup, chip=chip, init=init, clobber=clobber, $
  flat=flat, arc=arc, slitflat=slitflat, proc=proc, dotrace=dotrace, $
  skysub=skysub, extract=extract, combine=combine, calibrate=calibrate, $
  dostandards=dostandards, makesens=makesens, final_spec1d=final_spec1d, $
  tarballs=tarballs
; jm09sepucsd - uber-wrapper to reduce *all* the hizea/hires data 
    
;   hires_rslvall
    
    hizeapath = getenv('HIZEA_HIRES_DIR')+'/'
    spec1dpath = hizeapath+'spec1d/'
    if (n_elements(night) ne 0L) then $
      allnight = night else $
;     allnight = '09dec11'
      allnight = ['09dec11','09dec10','09aug13','09aug12']
    nnight = n_elements(allnight)

    xtoler = 0.002
    chip = [1,2] ; [1=blue, 2=green, 3=red]
    
; ##################################################
; data reduction stuff    
    for inight = 0, nnight-1 do begin
       datapath = hizeapath+allnight[inight]+'/'
       print & splog, '##################################################'
       splog, datapath
; define the various setups, which correspond to various
; cross-disperser angles (see hires_summ.txt)
       case allnight[inight] of
          '09dec11': allsetup = [1,2,3,4]
          '09dec10': allsetup = [1,2,3]
          '09aug13': allsetup = [1,2]
          '09aug12': allsetup = [1,2]
          else: stop
       endcase
       nsetup = n_elements(allsetup)
; ##################################################
; initialize all the directories and run the preprocessing required by
; the HIRES pipeline; this only should be run once!!
       if keyword_set(init) then begin
; now run the auto-id routine
          if file_test(datapath+'rawhires.fits',/regular) and $
            (keyword_set(clobber) eq 0) then begin
             rawhires = hires_ar(datapath+'rawhires.fits') 
          endif else begin
             pushd, datapath
             hires_strct, rawhires, outfil=datapath+'rawhires.fits'
             popd
          endelse 
; each night needs a customized routine to edit the initial HIRES
; structure; for example, to fix headers, reject crappy spectra, etc.
          case allnight[inight] of
             '09dec11': hires = edit_hires_09dec11(rawhires)
             '09dec10': hires = edit_hires_09dec10(rawhires)
             '09aug13': hires = edit_hires_09aug13(rawhires)
             '09aug12': hires = edit_hires_09aug12(rawhires)
             else: stop
          endcase
; run some final setup steps and then write out 
          hires_setup, hires, xtoler=xtoler, outfil=datapath+'hires_summ.txt'
; calculate the gain, for each setup and then write out
;         for jj = 0, nsetup-1 do hires_findgain, hires, allsetup[jj];, chip
          hires_wrstrct, hires, fits=datapath+'hires.fits', $
            outfil=datapath+'hires.list'
          continue
       endif 
; ##################################################
; now reduce each night, for each setup
       if (keyword_set(final_spec1d) eq 0) and (keyword_set(tarballs) eq 0) then begin
          for jj = 0, nsetup-1 do begin
             hizea_reduce_night, datapath, allsetup[jj], chip=chip, clobber=clobber, $
               flat=flat, arc=arc, slitflat=slitflat, proc=proc, dotrace=dotrace, $
               skysub=skysub, extract=extract, combine=combine, calibrate=calibrate, $
               dostandards=dostandards, makesens=makesens
          endfor 
       endif
; ##################################################
; build the final 1D spectra
       if keyword_set(final_spec1d) then begin
          hires = hires_ar(datapath+'hires.fits')
          objindx = where(hires.type eq 'OBJ',nobj)

          objinfo = struct_trimtags(hires[objindx],select=['obj_id',$
            'obj','frame','chip','xdangl','ra','dec','exp','am','img_root'])
          info = struct_addtags(replicate({galaxy: '', z: 0.0, spec1dfile: ''},nobj),objinfo)

          for iobj = 0, nobj-1 do begin
; August 2009
             if strmatch(info[iobj].obj,'*1558*') then begin
                info[iobj].galaxy = 'J1558+3957'
                info[iobj].z = 0.402242
             endif
             if strmatch(info[iobj].obj,'*1658*') then begin
                info[iobj].galaxy = 'J1658+2354'
                info[iobj].z = 0.498282
             endif
             if strmatch(info[iobj].obj,'*2116*') then begin
                info[iobj].galaxy = 'J2116-0634'
                info[iobj].z = 0.728
             endif
             if strmatch(info[iobj].obj,'*2140*') then begin
                info[iobj].galaxy = 'J2140+1209'
                info[iobj].z = 0.75097
             endif
             if strmatch(info[iobj].obj,'*2256*') then begin
                info[iobj].galaxy = 'J2256+1542'
                info[iobj].z = 0.727
             endif
; December 2009
             if strmatch(info[iobj].obj,'*0826*') then begin
                info[iobj].galaxy = 'J0826+4305'
                info[iobj].z = 0.603
             endif
             if strmatch(info[iobj].obj,'*0901*') then begin
                info[iobj].galaxy = 'J0901+0314'
                info[iobj].z = 0.4586
             endif
             if strmatch(info[iobj].obj,'*0905*') then begin
                info[iobj].galaxy = 'J0905+5759'
                info[iobj].z = 0.7114
             endif
             if strmatch(info[iobj].obj,'*0944*') then begin
                info[iobj].galaxy = 'J0944+0930'
                info[iobj].z = 0.5138
             endif
          endfor 

; coadd all the spectra and make a QAplot          
          qafile = datapath+'spec1d/qa_'+allnight[inight]+'_spec1d.ps'
          hizea_coadd_spec1d, hires[objindx], info, $
            datapath=datapath, qafile=qafile
          qaplot_hizea_spec1d, datapath=datapath, qafile=qafile
       endif 
    endfor 

; ##################################################
; build the final 1D spectra
;   if keyword_set(final_spec1d) then begin
;      hizea_final_spec1d, hizeapath+allnight+'/spec1d/', spec1dpath, fluxed=0
;      hizea_final_spec1d, hizeapath+allnight+'/spec1d/', spec1dpath, fluxed=1
;   endif
       
; ##################################################
; make tarballs
    if keyword_set(tarballs) then begin
; 1d spectra from the individual nights
       for ii = 0L, nnight-1L do begin
          tarname = hizeapath+'tar/'+allnight[ii]+'_spec1d.tar.gz'
          pushd, hizeapath+allnight[ii]+'/spec1d/'
          flist = [file_search('*.fits.gz'),file_search('*.ps.gz')]
          spawn, 'tar czvf '+tarname+' '+strjoin(flist,' '), /sh
          popd
       endfor
;; final coadded spectra
;       tarname = hizeapath+'tar/final_spec1d.tar.gz'
;       pushd, hizeapath+'spec1d/'
;       flist = [file_search('*.fits'),file_search('*.ps.gz')]
;       spawn, 'tar czvf '+tarname+' '+strjoin(flist,' '), /sh
;       popd
    endif
       
return
end

;;; stupidly hard-coded for now!!
;;          allfiles = file_search(datapath+'FSpec/*.fits.gz',count=nfiles)
;;          pad = 30
;;          npix = 2000-pad*2
;;          for ifile = 0, nfiles-1 do begin
;;             out1 = {file: '', object: '', ordr: 0, keep: 1, $
;;               wave: fltarr(npix), flux: fltarr(npix)}
;;             splog, 'Reading '+allfiles[ifile]
;;             fspec = mrdfits(allfiles[ifile],1,/silent)
;;             these = where(fspec.phys_ordr gt 0,nthese)
;;
;;             out = replicate(out1,nthese)
;;             out.file = strtrim(allfiles[ifile],2)
;;             out.object = strtrim(fspec.field,2)
;;             out.ordr = fspec.phys_ordr[these]
;;             for jj = 0, nthese-1 do begin
;;                good = where(fspec.wave[*,these[jj]] gt 0.0,ngood)
;;;               if fspec.phys_ordr[these[jj]] eq 91 then stop
;;;               good = good[where(good le npix-1)]
;;;               ngood = n_elements(good)
;;                if (ngood-2*pad lt npix) then out[jj].keep = 0 else begin
;;                   out[jj].wave = fspec.wave[good[pad:ngood-pad-1],these[jj]]
;;                   out[jj].flux = fspec.fx[good[pad:ngood-pad-1],these[jj]]
;;                   crap = where(out[jj].wave eq 0.0)
;;                   if (crap[0] ne -1) then stop
;;                endelse
;;             endfor
;;             if (ifile eq 0) then allout = out else allout = [allout,out]
;;          endfor 
;;; J1558+3957
;;          these = where(strmatch(allout.object,'*J1558*'))
;;          outfile = datapath+'spec1d/J1558+3957_hires.fits'
;;          final = allout[these]
;;          final = final[reverse(sort(final.ordr))]
;;          final = final[where(final.keep)]
;;          im_mwrfits, final, outfile
;;; J2140+1209
;;          these = where(strmatch(allout.object,'*J2140*'))
;;          outfile = datapath+'spec1d/J2140+1209_hires.fits'
;;          final = allout[these]
;;          final = final[reverse(sort(final.ordr))]
;;          final = final[where(final.keep)]
;;          im_mwrfits, final, outfile
;;;; sort by unique object number
;;;          allobj = strtrim(allout.object,2)
;;;          obj = allobj[uniq(allobj,sort(allobj))]
;;;          nobj = n_elements(obj)
;;;          for kk = 0, nobj-1 do begin
;;;             these = where(obj[kk] eq allobj)
;;;             final = allout[these]
;;;             sort = sort(final.ordr)
;;;          endfor
