function j1506_path, isedfit=isedfit, montegrids=montegrids
; jm12mar19ucsd

    path = getenv('HIZEA_DATA')+'/j1506/'
    if keyword_set(isedfit) then path = path+'isedfit/'
    if keyword_set(montegrids) then path = path+'montegrids/'
    
return, path
end
