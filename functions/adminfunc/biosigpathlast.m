function biosigpathlast

str2doublepath = fileparts( which('str2double') );
sopenpath      = fileparts( which('sopen') );
if strcmp(str2doublepath,sopenpath)
    rmpath(str2doublepath);
    str2doublepat2 = fileparts( which('str2double') );
    addpath(str2doublepath,'-begin');
    addpath(str2doublepat2,'-begin');
end

