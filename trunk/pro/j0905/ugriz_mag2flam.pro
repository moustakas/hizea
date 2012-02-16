
function ugriz_mag2flam, mag, mag_err, flux_err=flux_err

if n_elements(mag) ne 5 then begin
  print, 'mag must be a 5 element array corresponding to ugriz'
  return, 0
endif

flux = fltarr(5)
flux[0] = sdss_mag2flam(mag[0], mag_err[0], 'u', flux_err=uerr)
flux[1] = sdss_mag2flam(mag[1], mag_err[1], 'g', flux_err=gerr)
flux[2] = sdss_mag2flam(mag[2], mag_err[2], 'r', flux_err=rerr)
flux[3] = sdss_mag2flam(mag[3], mag_err[3], 'i', flux_err=ierr)
flux[4] = sdss_mag2flam(mag[4], mag_err[4], 'z', flux_err=zerr)

flux_err = [uerr, gerr, rerr, ierr, zerr]

return, flux

end
