function biosigpathfirst

str2doublepath = fileparts( which('str2double') );
sopenpath      = fileparts( which('sopen') );
if ~strcmp(str2doublepath,sopenpath)
    rmpath(sopenpath);
    addpath(sopenpath,'-begin');
end


