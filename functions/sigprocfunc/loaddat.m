% loaddat() - loading neuroscan format data file into matlab.
%
% Usage:
%   >> [typeeeg, rt, response, n] = loaddat( filename );
%
% Inputs:
%   filename     - input Neuroscan .dat file
%
% Outputs:
% 	typeeeg	 - type of the sweeps
%	rt	 - reaction time of the subject
%	response - response of the subject
%	n	 - number of sweeps
%
% See also: loadeeg()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [tmptypeeeg, tmprt, tmpresponseeeg, n] = loaddat( FILENAME)

if nargin<1 
	fprintf('Not enought arguments\n'); 
	help loaddat 
	return;
end;

BOOL='int16';
ULONG='int32'; 
FLOAT='float32';

fid=fopen(FILENAME,'r','ieee-le');
if fid<0
	fprintf(2,['Error LOADEEG: File ' FILENAME ' not found\n']);  
	return;
end;

% skip the first 20 lines
% -----------------------
for index=1:20	fgetl(fid); end;

% read the matrix
% ---------------
tmpMat 		    = fscanf(fid, '%f', [5, inf]);
tmptypeeeg  	= tmpMat(3,:);
tmpresponseeeg  = tmpMat(4,:);
tmprt       	= tmpMat(5,:);
n        	    = size( tmpMat, 2);

fclose(fid); 
return;


















