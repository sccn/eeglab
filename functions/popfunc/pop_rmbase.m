% pop_rmbase() - remove baseline from an epoched dataset.
%
% Usage:
%   >> outeeg = pop_rmbase( eeg, timerange, pointrange);
%
% Inputs:
%   eeg        - input dataset
%   timerange  - timerange for baseline [min max] in milliseconds
%   pointrange - pointrange for baseline [min max]. Timerange must be
%                empty. Timerange is prioritary over pointrange.
%
% Outputs:
%   outeeg     - output dataset
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: rmbase(), eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.2  2002/08/12 02:29:51  arno
% inputdlg2
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 

function [EEG, com] = pop_rmbase( EEG, timerange, pointrange);

com ='';
if isempty(EEG.data)
    disp('Pop_rmbase error: cannot remove baseline of empty dataset'); return;
end;    
if nargin < 1
	help pop_rmbase;
	return;
end;	
if nargin < 2
	% popup window parameters
	% -----------------------
   promptstr    = {'Enter time limits (in ms) for baseline:',...
         		   strvcat('Enter point range for baseline (ex:1:frames):', '(overwrite previous option)') };
	inistr       = { [num2str(EEG.xmin*1000) ' 0'], '' };
	result       = inputdlg2( promptstr, 'Baseline removal -- pop_rmbase()', 1,  inistr, 'pop_rmbase');
	size_result  = size( result );
	if size_result(1) == 0 return; end;

	% decode parameters
	% -----------------
	if isempty( result{1} ) & isempty( result{2} ), return; end;
    timerange = eval( [ '[' result{1} ']' ] );
    pointrange = eval( [ '[' result{2} ']' ] );
end;

if ~isempty(timerange)
    pointrange = round((timerange(1)/1000-EEG.xmin)*EEG.srate+1):round((timerange(2)/1000-EEG.xmin)*EEG.srate);
end;	

if ~isempty(timerange) & (timerange(1) < EEG.xmin*1000) & (timerange(2) > EEG.xmax*1000)
   error('Wrong time range');
end;
if (min(pointrange) < 1) | (max( pointrange ) > EEG.pnts)  
   error('Wrong point range');
end;

fprintf('Removing baseline...\n');
EEG.data = rmbase( EEG.data(:,:), EEG.pnts, pointrange );
EEG.data = reshape( EEG.data, EEG.nbchan, EEG.pnts, EEG.trials);

if ~isempty(timerange)
	com = sprintf('%s = pop_rmbase( %s, [%s]);', inputname(1), inputname(1), ...
			num2str(timerange));
else
	com = sprintf('%s = pop_rmbase( %s, [], [%s]);', inputname(1), inputname(1), ...
			num2str(pointrange));
end;			

return;
