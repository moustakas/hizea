pro sigmasfr_hstmontage
; jm12apr05ucsd - make a montage
    
    path = sigmasfr_path()
    outfile = path+'hst_montage.png'

    hst = rsex(path+'hst_sample.dat')
    lir = mrdfits(path+'sigmasfr_lir.fits.gz',1)
    good = where(lir.sfr_chary gt -900)
    hst = hst[good]
    lir = lir[good]
    
    sigmasfr = alog10(lir.sfr_chary/(!pi*hst.r_e^2)/2.0)

; sort by sigmasfr and then vout    
    srt = multisort(reverse(sigmasfr),reverse(hst.vout))
    sigmasfr = sigmasfr[srt] & hst = hst[srt]
    niceprint, hst.galaxy, sigmasfr, hst.vout

    pushd, path+'jpeg'
    alljpeg = file_search(strlowcase(strtrim(hst.galaxy,2))+'.png')
    cmd = 'montage -bordercolor black -borderwidth 1 '+ $
      '-tile 5x5 -geometry +0+0 -quality 100 -resize 200x200 '+$
      strjoin(alljpeg,' ')+' '+outfile
    spawn, cmd, /sh
    popd

return
end
    
