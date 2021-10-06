function biosigpathfirst
% Add BIOSIG at the beginning of the path 

str2doublepath = fileparts( which('str2double') );
sopenpath      = fileparts( which('sopen') );
if ~strcmp(str2doublepath,sopenpath)
    addpath(sopenpath,'-begin');
end


