% pop_rmbase() - remove channel baseline means from an epoched or 
%                continuous EEG dataset. Calls rmbase().
% Usage:
%   >> OUTEEG = pop_rmbase( EEG ); % pop up an interactive arg entry window
%   >> OUTEEG = pop_rmbase( EEG, timerange, pointrange); % call rmbase()
%
% Graphic interface:
%    "Baseline latency range" - [edit box] Latency range for the baseline in ms.
%                               Collects the 'timerange' command line input.
%                               Empty or [] input -> Use whole epoch as baseline
%    "Baseline points vector" - [edit box] Collects the 'pointrange' command line 
%                               option (below). (Overwritten by 'timerange'). 
%                               Empty or [] input -> Use whole epoch as baseline
% Inputs:
%   EEG        - Input dataset
%   timerange  - [min_ms max_ms] Baseline latency range in milliseconds.
%                                Empty or [] input -> Use whole epoch as baseline
%   pointrange - [min:max]       Baseline points vector (overwritten by timerange).
%                                Empty or [] input -> Use whole epoch as baseline
% Outputs:
%   OUTEEG     - Output dataset
%
% Note: If dataset is continuous, channel means are removed separately 
%       for each continuous data region, respecting 'boundary' events 
%       marking boundaries of excised or concatenated data portions.
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
% Revision 1.16  2005/02/27 02:59:22  scott
% make entries in ms as per documentation
% added ability to use empty inputs -> use whole epoch for baseline
%
% Revision 1.15  2004/12/16 23:21:13  arno
% header in ms
%
% Revision 1.14  2004/12/16 22:59:32  hilit
% fixing some typo with 'pointrange' and separating the | condition in the if
% statement
%
% Revision 1.13  2004/08/18 16:45:01  arno
% ms -> s
%
% Revision 1.12  2003/07/18 14:50:28  scott
% commenting, editting help message
%
% Revision 1.11  2003/05/12 15:39:51  arno
% debug command line call
%
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
    promptstr    = {'Baseline latency range (min_ms max_ms) ([] = whole epoch):',...
         		   strvcat('Else, baseline points vector (ex:1:56) ([] = whole epoch):', ...
                           '(overwritten by latency range above).') };
	inistr     = { [num2str(EEG.xmin*1000) ' 0'], '' }; % latency range in ms
	result       = inputdlg2( promptstr, 'Epoch baseline removal -- pop_rmbase()', ...
                                  1,  inistr, 'pop_rmbase');
	size_result  = size( result );
	if size_result(1) == 0 return; end;

	% decode parameters
	% -----------------
	if numel(result) < 2 | ((isempty(result{1}) | strcmp(result{1},'[]') ) ...
                           & (isempty(result{2}) | strcmp(result{2},'[]')))
       timerange = [num2str(EEG.xmin*1000) num2str(EEG.xmax*1000)]; % whole epoch latency range
       fprintf('pop_rmbase(): using whole epoch as baseline.\n');

       % fprintf('pop_rmbase(): baseline limits must be specified.\n');
       % return; end;
    else
      timerange = eval( [ '[' result{1} ']' ] );
      pointrange = eval( [ '[' result{2} ']' ] );
    end
elseif nargin < 2 & EEG.trials == 1
	% popup window parameters
	% -----------------------
    resp = questdlg2(strvcat('Remove baseline of each channel (in', ...
             'each unbroken portion of continuous data)'), 'pop_rmbase', 'Cancel', 'Ok', 'Ok');
    if strcmpi(resp, 'Cancel'), return; end;
    timerange = [];
    pointrange = [1:EEG.pnts];
end;

if exist('pointrange') ~= 1
    if ~isempty(timerange) & (timerange(1) < EEG.xmin*1000) & (timerange(2) > EEG.xmax*1000)
        error('pop_rmbase(): Bad time range');
    end;
    pointrange = round((timerange(1)/1000-EEG.xmin)*EEG.srate+1):round((timerange(2)/1000-EEG.xmin)*EEG.srate);
end;

if isempty(pointrange)
    if ~isempty(timerange) & (timerange(1) < EEG.xmin*1000) & (timerange(2) > EEG.xmax*1000)
        error('pop_rmbase(): Bad time range');
    end;
    pointrange = round((timerange(1)/1000-EEG.xmin)*EEG.srate+1):round((timerange(2)/1000-EEG.xmin)*EEG.srate);
end;	

if (min(pointrange) < 1) | (max( pointrange ) > EEG.pnts)  
   error('pop_rmbase(): Wrong point range');
end;

fprintf('pop_rmbase(): Removing baseline...\n');
%
% Respect excised data boundaries if continuous data
% ---------------------------------------------------
if EEG.trials == 1 & ~isempty(EEG.event) ...
                     & isfield(EEG.event, 'type') ...
                        & isstr(EEG.event(1).type)
	boundaries = strmatch('boundary', {EEG.event.type});
	if ~isempty(boundaries)
        fprintf('Pop_rmbase(): finding continuous data discontinuities\n');
        boundaries = cell2mat({EEG.event(boundaries).latency})-0.5-pointrange(1)+1;
        boundaries(find(boundaries>=pointrange(end)-pointrange(1))) = [];
        boundaries(find(boundaries<1)) = [];
        boundaries = [0 boundaries pointrange(end)-pointrange(1)];
        for index=1:length(boundaries)-1
            tmprange = [boundaries(index)+1:boundaries(index+1)];
            EEG.data(:,tmprange) = rmbase( EEG.data(:,tmprange), length(tmprange), ...
                                               [1:length(tmprange)]);
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
