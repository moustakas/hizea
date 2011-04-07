pro write_hizea_casjobs_input
; jm10dec22ucsd - write the coordinate file that needs to be uploaded
;   to the MAST/casjobs database; see
;   /Users/ioannis/home/research/projects/ages/catalogs/galex/README
;   for more details

    iopath = getenv('HIZEA_DATA')+'/'
    parent = mrdfits(iopath+'sdss/hizea_simplefit_dr7.fit',1)
    ngal = n_elements(parent)

    out = struct_addtags(replicate({id: 0L},ngal),$
      struct_trimtags(parent,select=['ra','dec']))
    out.id = lindgen(ngal)
    
    outfile = iopath+'galex/hizea_casjobs_input.dat'
    openw, lun, outfile, /get_lun
    printf, lun, '# id ra dec'
    struct_print, out, lun=lun, ddigit=12, /no_head
    free_lun, lun

return
end
    
