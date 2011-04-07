function casjobs_remove_double, cat
; most of the fields in the casjobs catalog are unnecessarily double
; precision; change them to float
    tags = strlowcase(tag_names(cat))
    ntags = n_elements(tags)
    for ii = 0, ntags-1 do begin
       type = size(cat[0].(ii),/type)
       if (type eq 5) then begin
          if (tags[ii] eq 'alpha_j2000') or (tags[ii] eq 'delta_j2000') or $
            (tags[ii] eq 'ra') or (tags[ii] eq 'dec') then type = 5 else type = 4
       endif
       if (type eq 3) then begin
          if (tags[ii] eq 'hizea_id') or (tags[ii] eq 'tile') or $
            (tags[ii] eq 'objid') then type = 3 else type = 2
       endif
       arr = (make_array(1,type=type))[0]
;      print, tags[ii], type, arr
       if (ii eq 0) then outcat = create_struct(tags[ii],arr) else $
         outcat = create_struct(outcat,tags[ii],arr)
    endfor
    outcat = replicate(outcat,n_elements(cat))
    for ii = 0, ntags-1 do outcat.(ii) = cat.(ii)
return, outcat
end

pro build_hizea_galex_catalog, out_galex, gr=gr, clobber=clobber
; jm10apr30ucsd - build a line-matched GALEX catalog for the HIZEA
; sample using the CasJobs output (see the README in the
; $HIZEA_REDUX/galex directory)

    if (n_elements(gr) eq 0) then gr = 'gr6'
    hizeapath = getenv('HIZEA_DATA')+'/sdss/'

; see WRITE_HIZEA_CASJOBS_INPUT for how the input catalog was written
    hizea = mrdfits(hizeapath+'hizea_simplefit_dr7.fit',1)
    ngal = n_elements(hizea)
    id = lindgen(ngal)

; output filename    
    outfile = hizeapath+'hizea_galex_'+gr+'.fits'
    if file_test(outfile+'.gz') and (keyword_set(clobber) eq 0) then begin
       splog, 'Output file '+outfile+' exists; use /CLOBBER'
       return
    endif
    
; read the casjobs output; for some reason some input/output objects
; are repeated (not sure why), so remove them here
    incat = mrdfits(hizeapath+'hizea_galex_'+gr+'_casjobs.fits.gz',1)
    incat = casjobs_remove_double(incat)
    incat = struct_addtags(replicate({galex_object_position: -999L},$
      n_elements(incat)),temporary(incat))
    incat.galex_object_position = lindgen(n_elements(incat))
    
    allid = strtrim(incat.id,2)+strtrim(incat.objid,2)
    keep = uniq(allid,sort(allid))
    splog, 'Removing '+strtrim(n_elements(incat)-n_elements(keep),2)+$
      '/'+strtrim(n_elements(incat),2)+' repeated casjobs entries!'
    incat = incat[keep]

    outcat1 = galex_resolve_casjobs_duplicates(incat,idtag='id')
    best = where(outcat1.isbest,nbest)
    outcat = outcat1[best]
    splog, 'Found '+strtrim(nbest,2)+'/'+strtrim(ngal,2)+$
      ' galaxies with good GALEX detections!'

    match, id, outcat.id, m1, m2
    srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
    out_galex = im_empty_structure(outcat[0],$
      empty_value=-999.0,ncopies=ngal)
    out_galex[m1] = outcat[m2]

; add ID, RA, DEC from the original HIZEA structure to the output
; catalog
    out_galex.id = id
    out_galex.ra = hizea.ra
    out_galex.dec = hizea.dec
    im_mwrfits, out_galex, outfile, /clobber

return
end
    
