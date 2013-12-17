pro build_hizea_chandra_sample, clobber=clobber
; jm13aug31siena - gather everything we need for the chandra sample

    chandrapath = getenv('HIZEA_DATA')+'/chandra/'
    mmtpath = chandrapath+'mmt/'

    allspec = file_search(mmtpath+'*.fits',count=nspec)
    mmt = mrdfits(allspec[0],1,/silent,columns=['short_name','z'])
    for ii = 1, nspec-1 do mmt = [mmt,mrdfits(allspec[ii],1,$
      /silent,columns=['short_name','z'])]
    
;    gal = [$
;      'J0826+4305',$
;      'J0944+0930',$
;      'J1104+5946',$
;      'J1359+5137',$
;      'J1506+5402',$
;      'J1506+6131',$
;      'J1558+3957',$
;      'J1613+2834',$
;      'J1634+4619',$
;      'J1713+2817',$
;      'J2118+0017',$
;      'J2140+1209']
;;   id = lindgen(n_elements(gal))+1

; get the photometry
    cat = mrdfits(sigmasfr_path()+'sigmasfr_photometry.fits.gz',1)
;   cat = struct_addtags(replicate({id: 0},n_elements(cat)),cat)

    match, strtrim(mmt.short_name,2), strtrim(cat.galaxy,2), m1, m2
    srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
    outcat = cat[m2]
;   outcat.id = id[m1]
    niceprint, outcat.galaxy, mmt[m1].short_name, outcat.z, mmt[m1].z
    
; sort by redshift!
    outcat = outcat[sort(outcat.z)]
    
    im_mwrfits, outcat, chandrapath+'hizea_chandra_photometry.fits', $
      clobber=clobber
    
return
end
