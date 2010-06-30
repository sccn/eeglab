% helpforexe() - Write help files for exe version
%
% Usage:
%   histtoexe( mfile, folder)
%
% Inputs:
%   mfile  - [cell of string] Matlab files with help message
%   folder - [string] Output folder
%
% Output:
%   text files name help_"mfile".m
%
% Author: Arnaud Delorme, 2006
%
% See also: eeglab() 

% Copyright (C) 2006 Arnaud Delorme
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

function helpforexe( funct, fold );

if nargin <1
	help histforexe;
	return;
end;
nonmatlab = 0;

% write all files
% ---------------
for index = 1:length(funct)
    doc1 = readfunc(funct{index}, nonmatlab);
    
    fid = fopen( fullfile(fold, [ 'help_' funct{index} ]), 'w');
    for ind2 = 1:length(doc1)
        if isempty(doc1{ind2}) fprintf(fid, [ 'disp('' '');\n' ]);
        else                   fprintf(fid, [ 'disp(' vararg2str({ doc1{ind2} }) ');\n' ]);
        end;
    end;
    fclose(fid);
    %fprintf(fid, 'for tmpind = 1:length(tmptxt), if isempty(tmptxt{tmpind}), disp('' ''); else disp(tmptxt{tmpind}); end; end; clear tmpind tmptxt;\n' );
end;

% try all functions
% -----------------
tmppath = pwd;
cd(fold);
for index = 1:length(funct)
    evalc([ 'help_' funct{index}(1:end-2) ]);
end;
cd(tmppath);
return;

%-------------------------------------
function [doc] = readfunc(funct, nonmatlab)

doc = {};
if nonmatlab	
	fid = fopen( funct, 'r');
else
	if findstr( funct, '.m')
		fid = fopen( funct, 'r');
	else
		fid = fopen( [funct '.m'], 'r');
	end;
end;

if fid == -1
	error('File not found');
end;

sub = 1;
try, 
    if ~isunix, sub = 0; end;
catch, end;

if nonmatlab
	str = fgets( fid );
	while ~feof(fid)
		str = deblank(str(1:end-sub));
        doc = { doc{:} str(1:end) };    
        str = fgets( fid );
	end;
else
	str = fgets( fid );
	while (str(1) == '%')
		str = deblank(str(1:end-sub));
        doc = { doc{:} str(2:end) };    
		str = fgets( fid );
	end;
end;
fclose(fid);
