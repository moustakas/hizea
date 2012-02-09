function j0905_path, isedfit=isedfit, montegrids=montegrids
; jm12feb08ucsd

    hizeapath = getenv('HIZEA_DATA')+'/j0905/'

    if keyword_set(isedfit) then hizeapath = hizeapath+'isedfit/'
    if keyword_set(montegrids) then hizeapath = hizeapath+'montegrids/'
    
return, hizeapath
end
