function biosigpathlast
% Add BIOSIG at the end of the path to avoid overloading Matlab functions

str2doublepath  = fileparts( which('str2double') );
str2numpath     = fileparts( which('str2num') );
if ~strcmp(str2doublepath,str2numpath)
    addpath(str2numpath,'-begin');
end

