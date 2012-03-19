function sigmasfr_path, isedfit=isedfit, montegrids=montegrids
; jm12mar19ucsd

    path = getenv('HIZEA_DATA')+'/sigmasfr/'
    if keyword_set(isedfit) then path = path+'isedfit/'
    if keyword_set(montegrids) then path = path+'montegrids/'
    
return, path
end
