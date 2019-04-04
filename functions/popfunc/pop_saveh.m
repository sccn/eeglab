% pop_saveh() - save the EEGLAB session command history stored in ALLCOM
%               or in the 'history' field of the current dataset
%
% Usage:
%   >> pop_saveh( ALLCOM, filename, filepath);
%   >> pop_saveh( EEG.history, filename, filepath);
%
% Inputs:
%   ALLCOM      - cell array of strings containing the EEGLAB command history 
%   EEG.history - history field of the current dataset
%   filename    - name of the file to save to (optional, default "eeglabhist.m"
%   filepath    - path of the file to save to (optional, default pwd)
%
% Author: Arnaud Delorme, CNL / Salk Institute, 22 March 2002
%
% See also: eegh(), eeglab()

% Copyright (C) 22 March 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% 01-25-02 reformated help & license -ad 
% 03-29-02 added update menu -ad 

function com = pop_saveh( allcoms, curfilename, curfilepath);

com = '';
if nargin < 1
	help pop_saveh;
	return;
end
	
if nargin < 3
	[curfilename, curfilepath] = uiputfile('eeglabhist.m', 'Save the EEGLAB session command history with .m extension -- pop_saveh()');
    drawnow;
	if curfilename == 0 return; end
end


fid = fopen( [ curfilepath '/' curfilename ], 'w');
if fid == -1
    error('pop_saveh(): Cannot open named file');
end;    
fprintf(fid, '%% EEGLAB history file generated on the %s\n', date);
fprintf(fid, '%% ------------------------------------------------\n');
if iscell(allcoms)
    disp('Saving the EEGLAB session command history...');
    for index = length(allcoms):-1:1
        fprintf(fid, '%s\n', allcoms{index});
    end
    fprintf(fid, 'eeglab redraw;\n');
else
    disp('Saving the current EEG dataset command history...');
    fprintf(fid, '%s\n', allcoms);
end
fclose(fid);

if iscell(allcoms) && nargout == 1
    com = sprintf('pop_saveh( %s, ''%s'', ''%s'');', inputname(1), curfilename, curfilepath);
else
    com = sprintf('pop_saveh( EEG.history, ''%s'', ''%s'');', curfilename, curfilepath);
end

return;
