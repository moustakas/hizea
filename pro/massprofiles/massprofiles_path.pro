function massprofiles_path, isedfit=isedfit, montegrids=montegrids, hst=hst
; jm12mar19ucsd

    path = getenv('HIZEA_DATA')+'/massprofiles/'
    if keyword_set(isedfit) then path = path+'isedfit/'
    if keyword_set(montegrids) then path = path+'montegrids/'
    if keyword_set(hst) then path = path+'hst/'
    
return, path
end
