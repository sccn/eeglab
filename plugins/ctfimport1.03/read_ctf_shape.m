function [shape] = read_ctf_shape(filename);

% READ_CTF_SHAPE reads headshape points and header information
% from a CTF *.shape teh accompanying *.shape_info file.
%
% Use as
%   [shape] = read_ctf_shape(filename)
% where filename should have the .shape extension 

% Copyright (C) 2003, Robert Oostenveld
%
% $Log: not supported by cvs2svn $
% Revision 1.1  2005/12/06 06:24:23  psdlw
% Alternative functions from the FieldTrip package, which is now released under GPL (so I assume these functions can be committed to the sourceforge cvs)
%
% Revision 1.1  2003/10/08 15:28:03  roberto
% *** empty log message ***
%

shape = read_ctf_ascii([filename '_info']);

if ~strcmp(shape.MRI_Info.COORDINATES, 'HEAD')
  warning('points on head shape are NOT in headcoordinates')
end

fid = fopen(filename, 'rb');
num = fscanf(fid, '%d', 1);
shape.pnt = fscanf(fid, '%f', inf);
shape.pnt = reshape(shape.pnt, [3 num])';
fclose(fid);

