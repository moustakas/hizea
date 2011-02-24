;pro irac_wrapper
;+ 
;  amd101208 -- procedure that runs irac_phot.pro
;
;-

  !p.thick=2
  !x.thick=2
  !y.thick=2
  !p.charthick=2

  yo = file_search('r*/ch*/pbcd/*maic.fits')

  flux = fltarr(n_elements(yo))

  close, /all
  openw, 1, 'irac_photo_prelim.txt'

  dfpsplot, 'tmp.ps', /color

    for i=0L, n_elements(yo)-1L do begin
      flux[i] = irac_phot(yo[i], /nocntrd)
      printf, 1, yo[i], flux[i], format='(a13,2x,f7.1)'
    endfor
 

  dfpsclose

  close, /all


  !p.thick=1
  !x.thick=1
  !y.thick=1
  !p.charthick=1

stop

end
