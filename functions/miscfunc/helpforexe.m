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
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function helpforexe( funct, fold );

if nargin <1
	help histforexe;
	return;
end
nonmatlab = 0;

% write all files
% ---------------
for index = 1:length(funct)
    doc1 = readfunc(funct{index}, nonmatlab);
    
    fid = fopen( fullfile(fold, [ 'help_' funct{index} ]), 'w');
    for ind2 = 1:length(doc1)
        if isempty(doc1{ind2}) fprintf(fid, [ 'disp('' '');\n' ]);
        else                   fprintf(fid, [ 'disp(' vararg2str({ doc1{ind2} }) ');\n' ]);
        end
    end
    fclose(fid);
    %fprintf(fid, 'for tmpind = 1:length(tmptxt), if isempty(tmptxt{tmpind}), disp('' ''); else disp(tmptxt{tmpind}); end; end; clear tmpind tmptxt;\n' );
end

% try all functions
% -----------------
tmppath = pwd;
cd(fold);
for index = 1:length(funct)
    evalc([ 'help_' funct{index}(1:end-2) ]);
end
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
	end
end

if fid == -1
	error('File not found');
end

sub = 1;
try, 
    if ~isunix, sub = 0; end
catch, end

if nonmatlab
	str = fgets( fid );
	while ~feof(fid)
		str = deblank(str(1:end-sub));
        doc = { doc{:} str(1:end) };    
        str = fgets( fid );
	end
else
	str = fgets( fid );
	while (str(1) == '%')
		str = deblank(str(1:end-sub));
        doc = { doc{:} str(2:end) };    
		str = fgets( fid );
	end
end
fclose(fid);
