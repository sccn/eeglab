function [tmptypeeeg, tmprt, tmpresponseeeg, n] = loaddat( FILENAME)

% eeg_load_scan_dat - Read Neuroscan .dat format
%
% usage  [typeeeg, rt, response, n] = loaddat( FILENAME)
%
% Inputs:
%   FILENAME     input Neuroscan .dat file
% Outputs:
% 	typeeeg		type of the sweeps
%	rt		reaction time of the subject
%	response	response of the subject
%	n		number of sweeps
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:52 $

% Licence:  GNU GPL, no implied or express warranty
% History:  05/2001, arno_delorme@salk.edu
%           - Version 1.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
tmprt       	= tmpMat(5,:) * 1000;
n        	    = size( tmpMat, 2);

fclose(fid); 
return;
