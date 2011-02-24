function hizea_spitzer_forage, hdr
    info = {$
      object:   sxpar(hdr,'OBJECT'),$
      gain:     sxpar(hdr,'GAIN'),$
      fluxconv: sxpar(hdr,'FLUXCONV'),$
      exptime:  sxpar(hdr,'EXPTIME'),$
      aorkey:   sxpar(hdr,'AORKEY')}
return, info
end

pro unpack_hizea_spitzer, imagelist
; jm10dec08ucsd - unpack the raw Spitzer mosaics, rename them, and
; organize them in a more logical way

    rawpath = getenv('HIZEA_DATA')+'/spitzer/rawdata/'
    mosaicpath = getenv('HIZEA_DATA')+'/spitzer/mosaics/'
    aor = file_basename(file_search(rawpath+'r*',count=naor))
    
    ch1list = file_search(rawpath+aor+'/ch1/pbcd/SPITZER_I1_'+$
      strmid(aor,1)+'_????_?_E???????_maic.fits',count=nch1)
    ch2list = file_search(rawpath+aor+'/ch2/pbcd/SPITZER_I2_'+$
      strmid(aor,1)+'_????_?_E???????_maic.fits',count=nch2)

    for ii = 0, naor-1 do begin
       ch1 = mrdfits(ch1list[ii],0,hdr1)
       ch2 = mrdfits(ch2list[ii],0,hdr2)
       
       ch1info = hizea_spitzer_forage(hdr1)
       ch2info = hizea_spitzer_forage(hdr2)

stop       
    endfor
    
    
stop
    
    
return
end
