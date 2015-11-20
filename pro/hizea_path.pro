function hizea_path, phot=phot, isedfit=isedfit, montegrids=montegrids, code=code
; jm11apr06ucsd - 
; jm15nov17siena - updated to the latest directory structure

    datapath = getenv('HIZEA_DATA')+'/'

    path = datapath
    if keyword_set(phot) then path = datapath+'phot/'
    if keyword_set(isedfit) then path = datapath+'isedfit/'
    if keyword_set(montegrids) then path = datapath+'isedfit/montegrids/'
    if keyword_set(code) then path = getenv('HIZEA_DIR')+'/'
    
return, path
end
