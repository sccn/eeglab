% pop_rmbase() - remove baseline from an epoched or continuous
%                eeg dataset.
%
% Usage:
%   >> OUTEEG = pop_rmbase( EEG ); % pop up interactive window
%   >> OUTEEG = pop_rmbase( EEG, timerange, pointrange);
%
% Graphical interface:
%    "Baseline time range" - [edit box] Same as the 'timerange' command
%                line input.
%    "Baseline points vector" - [edit box] Overwritten by the time 
%                limits option above. Same as the 'pointrange' command 
%                line input.
%
% Inputs:
%   EEG        - Input dataset
%   timerange  - Baseline time range [min_ms max_ms]. 
%   pointrange - Baseline points vector [min:max]. 
%                (Overwritten by time range).
%
% Outputs:
%   OUTEGG     - Output dataset
%
% Note: in the case of a continuous dataset, the baseline is removed
%       seperatelly for the non-continuous regions (for instance if
%       data was removed).
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
% Revision 1.10  2003/04/25 18:35:35  arno
% now point range overwrite time limits
% /
%
% Revision 1.9  2003/02/17 02:56:16  arno
% reformating text for new functionality in help2html
%
% Revision 1.8  2003/02/16 23:53:12  arno
% typo last
%
% Revision 1.7  2003/02/16 23:52:49  arno
% adding gui info
%
% Revision 1.6  2002/11/18 02:47:23  arno
% adding continuous data baseline subtraction
%
% Revision 1.5  2002/10/11 14:48:20  arno
% recomputing ICA activations
%
% Revision 1.4  2002/08/17 20:26:22  scott
% menu and help msg
%
% Revision 1.3  2002/08/13 18:19:13  arno
% strvcat instead of 10
%
% Revision 1.2  2002/08/12 02:29:51  arno
% inputdlg2
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 

function [EEG, com] = pop_rmbase( EEG, timerange, pointrange);

com ='';
if nargin < 1
    help pop_rmbase;
    return;
end;
if isempty(EEG.data)
    disp('pop_rmbase(): cannot remove baseline of an empty dataset'); return;
end;    
if nargin < 1
	help pop_rmbase;
	return;
end;	
if nargin < 2 & EEG.trials > 1
	% popup window parameters
	% -----------------------
    promptstr    = {'Baseline time range (min_ms max_ms):',...
         		   strvcat('Baseline points vector (ex:1:56):', '(Overwrite time limits).') };
	inistr       = { [num2str(EEG.xmin*1000) ' 0'], '' };
	result       = inputdlg2( promptstr, 'Epoch baseline removal -- pop_rmbase()', 1,  inistr, 'pop_rmbase');
	size_result  = size( result );
	if size_result(1) == 0 return; end;

	% decode parameters
	% -----------------
	if isempty( result{1} ) & isempty( result{2} ), return; end;
    timerange = eval( [ '[' result{1} ']' ] );
    pointrange = eval( [ '[' result{2} ']' ] );
elseif nargin < 2 & EEG.trials == 1
	% popup window parameters
	% -----------------------
    resp = questdlg2(strvcat('Remove baseline for all data (separatelly for', ...
                     'each portion of continuous data)'), 'pop_rmbase', 'Cancel', 'Ok', 'Ok');
    if strcmpi(resp, 'Cancel'), return; end;
    timerange = [];
    pointrange = [1:EEG.pnts];
end;

if exist('pointarange') ~= 1 | isempty(pointrange)
    if ~isempty(timerange) & (timerange(1) < EEG.xmin*1000) & (timerange(2) > EEG.xmax*1000)
        error('pop_rembase(): Bad time range');
    end;
    pointrange = round((timerange(1)/1000-EEG.xmin)*EEG.srate+1):round((timerange(2)/1000-EEG.xmin)*EEG.srate);
end;	
if (min(pointrange) < 1) | (max( pointrange ) > EEG.pnts)  
   error('Wrong point range');
end;

fprintf('pop_rmbase(): Removing baseline...\n');
% add boundaries if continuous data
% ----------------------------------
if EEG.trials == 1 & ~isempty(EEG.event) & isfield(EEG.event, 'type') & isstr(EEG.event(1).type)
	boundaries = strmatch('boundary', {EEG.event.type});
	if ~isempty(boundaries)
        fprintf('Pop_rmbase: finding continuous data discontinuities\n');
        boundaries = cell2mat({EEG.event(boundaries).latency})-0.5-pointrange(1)+1;
        boundaries(find(boundaries>=pointrange(end)-pointrange(1))) = [];
        boundaries(find(boundaries<1)) = [];
        boundaries = [0 boundaries pointrange(end)-pointrange(1)];
        for index=1:length(boundaries)-1
            tmprange = [boundaries(index)+1:boundaries(index+1)];
            EEG.data(:,tmprange) = rmbase( EEG.data(:,tmprange), length(tmprange), [1:length(tmprange)]);
        end;
    else
        EEG.data = rmbase( EEG.data(:,:), EEG.pnts, pointrange );    
    end;		
else 
    EEG.data = rmbase( EEG.data(:,:), EEG.pnts, pointrange );
end;

EEG.data = reshape( EEG.data, EEG.nbchan, EEG.pnts, EEG.trials);
EEG.icaact = [];

if ~isempty(timerange)
	com = sprintf('%s = pop_rmbase( %s, [%s]);', inputname(1), inputname(1), ...
			num2str(timerange));
else
	com = sprintf('%s = pop_rmbase( %s, [], %s);', inputname(1), inputname(1), ...
			vararg2str({pointrange}));
end;			

return;
