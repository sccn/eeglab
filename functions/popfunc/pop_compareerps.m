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
end;
if exist('plottitle') ~= 1
	plottitle = '';
end;
if nargin < 2
   % which set to save
	% -----------------
	promptstr    = { 'List of datasets to compare (ex: 1 3 4):' ...
	                 'Channels subset to consider ([]=all):' ...
	                 'Plot title ([]=automatic):' };
	inistr       = { '1' '' '' };
	result       = inputdlg2( promptstr, 'Compare dataset ERPs -- pop_compareerps()', 1,  inistr, 'pop_compareerps');
	if length(result) == 0 return; end;
	setlist   	 = eval( [ '[' result{1} ']' ] );
	chansubset   = eval( [ '[' result{2} ']' ] );
	if isempty( chansubset ), chansubset = 1:ALLEEG(setlist(1)).nbchan; end;
	plottitle    = result{3};
    if isempty(plottitle)
        plottitle = [ 'Compare datasets number (blue=first; red=sec.)' int2str(setlist) ];
    end;
    figure;
end;
try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end;
tracing = [];
for setindex = setlist
	tracing  = [ tracing squeeze(mean(ALLEEG(setindex).data,3))];
end;

% save channel names
% ------------------
plottopo( tracing, ALLEEG(setlist(1)).chanlocs, ALLEEG(setlist(1)).pnts, [ALLEEG(setlist(1)).xmin ALLEEG(setlist(1)).xmax 0 0]*1000, plottitle, chansubset );

com = sprintf('figure; pop_compareerps( %s, [%s], [%s], ''%s'');', inputname(1), num2str(setlist), num2str(chansubset), plottitle);
return;
