function hizea_path, sdss=sdss, mass=mass, isedfit=isedfit, montegrids=montegrids
; jm11apr06ucsd - 

    hizeapath = getenv('HIZEA_DATA')+'/'

    if keyword_set(sdss) then hizeapath = hizeapath+'sdss/'
    if keyword_set(mass) then hizeapath = hizeapath+'mass/'
    if keyword_set(isedfit) then hizeapath = hizeapath+'mass/isedfit/'
    if keyword_set(montegrids) then hizeapath = hizeapath+'mass/montegrids/'
    
return, hizeapath
end
