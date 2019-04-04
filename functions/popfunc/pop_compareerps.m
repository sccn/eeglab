% pop_compareerps() - Compare the (ERP) averages of two datasets.
%
% Usage:
%       >> pop_compareerps( ALLEEG, datasetlist, chansubset, title);
% Inputs:
%   ALLEEG      - array of datasets
%   datasetlist - list of datasets
%   chansubset  - vector of channel subset
%   title       - plot title
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab(), plottopo()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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
% 03-11-02 added empty ALLEEG check -ad
% 03-18-02 added channel subset -ad
% 03-18-02 added title -ad & sm
function com = pop_compareerps( ALLEEG, setlist, chansubset, plottitle);

com = '';
if nargin < 1
   help pop_compareerps;
   return;
end;   
if isempty(ALLEEG)
    error('pop_compareerps: cannot process empty sets of data');
end
if exist('plottitle') ~= 1
	plottitle = '';
end
if nargin < 2
   % which set to save
	% -----------------
	promptstr    = { 'List of datasets to compare (ex: 1 3 4):' ...
	                 'Channels subset to consider ([]=all):' ...
	                 'Plot title ([]=automatic):' };
	inistr       = { '1' '' '' };
	result       = inputdlg2( promptstr, 'Compare dataset ERPs -- pop_compareerps()', 1,  inistr, 'pop_compareerps');
	if length(result) == 0 return; end
	setlist   	 = eval( [ '[' result{1} ']' ] );
	chansubset   = eval( [ '[' result{2} ']' ] );
	if isempty( chansubset ), chansubset = 1:ALLEEG(setlist(1)).nbchan; end
	plottitle    = result{3};
    if isempty(plottitle)
        plottitle = [ 'Compare datasets number (blue=first; red=sec.)' int2str(setlist) ];
    end
    figure;
end
try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end
tracing = [];
for setindex = setlist
	tracing  = [ tracing squeeze(mean(ALLEEG(setindex).data,3))];
end

% save channel names
% ------------------
plottopo( tracing, ALLEEG(setlist(1)).chanlocs, ALLEEG(setlist(1)).pnts, [ALLEEG(setlist(1)).xmin ALLEEG(setlist(1)).xmax 0 0]*1000, plottitle, chansubset );

if nargout == 1
    com = sprintf('figure; pop_compareerps( ALLEEG, [%s], [%s], ''%s'');', num2str(setlist), num2str(chansubset), plottitle);
end
return;
