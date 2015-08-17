function massprofiles_path, isedfit=isedfit, montegrids=montegrids, $
  hst=hst, code=code
; jm12mar19ucsd

    datapath = getenv('HIZEA_DATA')
    projectpath = getenv('HIZEA_PROJECT')

    path = projectpath+'/massprofiles/'
    if keyword_set(isedfit) then path = projectpath+'/massprofiles/isedfit/'
    if keyword_set(montegrids) then path = projectpath+'/massprofiles/isedfit/montegrids/'
    if keyword_set(hst) then path = datapath+'/hst/'
    if keyword_set(code) then path = getenv('HIZEA_DIR')+'/'
    
return, path
end
