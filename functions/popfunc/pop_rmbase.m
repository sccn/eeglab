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
%                               option (below). (Overwritten by 'timerange' above). 
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

function [EEG, com] = pop_rmbase( EEG, timerange, pointrange);

com ='';
if nargin < 1
    help pop_rmbase;
    return;
end;
if isempty(EEG(1).data)
    disp('pop_rmbase(): cannot remove baseline of an empty dataset'); return;
end;    
if nargin < 1
	help pop_rmbase;
	return;
end;	
if nargin < 2 & EEG(1).trials > 1
	% popup window parameters
	% -----------------------
    promptstr    = {'Baseline latency range (min_ms max_ms) ([] = whole epoch):',...
         		   strvcat('Else, baseline points vector (ex:1:56) ([] = whole epoch):', ...
                           '(overwritten by latency range above).') };
	inistr     = { [num2str(EEG(1).xmin*1000) ' 0'], '' }; % latency range in ms
	result       = inputdlg2( promptstr, 'Epoch baseline removal -- pop_rmbase()', ...
                                  1,  inistr, 'pop_rmbase');
	size_result  = size( result );
	if size_result(1) == 0 return; end;

	% decode parameters
	% -----------------
	if numel(result) < 2 | ((isempty(result{1}) | strcmp(result{1},'[]') ) ...
                           & (isempty(result{2}) | strcmp(result{2},'[]')))
       timerange = [num2str(EEG(1).xmin*1000) num2str(EEG(1).xmax*1000)]; % whole epoch latency range
       fprintf('pop_rmbase(): using whole epoch as baseline.\n');

       % fprintf('pop_rmbase(): baseline limits must be specified.\n');
       % return; end;
    else
      timerange = eval( [ '[' result{1} ']' ] );
      pointrange = eval( [ '[' result{2} ']' ] );
    end
elseif nargin < 2 & EEG(1).trials == 1
	% popup window parameters
	% -----------------------
    resp = questdlg2(strvcat('Remove mean of each data channel'), 'pop_rmbase', 'Cancel', 'Ok', 'Ok');
    if strcmpi(resp, 'Cancel'), return; end;
    timerange = [];
    pointrange = [1:EEG(1).pnts];
end;

% process multiple datasets
% -------------------------
if length(EEG) > 1
    [ EEG com ] = eeg_eval( 'pop_rmbase', EEG, 'warning', 'on', 'params', ...
                            { timerange pointrange } );
    return;
end;

if exist('pointrange') ~= 1 && ~isempty(timerange)
    if (timerange(1) < EEG.xmin*1000) & (timerange(2) > EEG.xmax*1000)
        error('pop_rmbase(): Bad time range');
    end;
    pointrange = round((timerange(1)/1000-EEG.xmin)*EEG.srate+1):round((timerange(2)/1000-EEG.xmin)*EEG.srate);
end;

if isempty(timerange)
    timerange = [ EEG(1).xmin*1000 EEG(1).xmax*1000];
end;

if exist('pointrange') ~= 1 || isempty(pointrange)
    if ~isempty(timerange) && (timerange(1) < EEG.xmin*1000) & (timerange(2) > EEG.xmax*1000)
        error('pop_rmbase(): Bad time range');
    end;
    pointrange = round((timerange(1)/1000-EEG.xmin)*EEG.srate+1):round((timerange(2)/1000-EEG.xmin)*EEG.srate);
end;	

if ~isempty(pointrange) && ((min(pointrange) < 1) || (max( pointrange ) > EEG.pnts))
   error('pop_rmbase(): Wrong point range');
end;

fprintf('pop_rmbase(): Removing baseline...\n');
%
% Respect excised data boundaries if continuous data
% ---------------------------------------------------
if EEG.trials == 1 && ~isempty(EEG.event) ...
                     && isfield(EEG.event, 'type') ...
                        && isstr(EEG.event(1).type)
    tmpevent = EEG.event;
	boundaries = strmatch('boundary', {tmpevent.type});
	if ~isempty(boundaries) % this is crashing
        fprintf('Pop_rmbase(): finding continuous data discontinuities\n');
        boundaries = round([ tmpevent(boundaries).latency ] -0.5-pointrange(1)+1);
        boundaries(boundaries>=pointrange(end)-pointrange(1)) = [];
        boundaries(boundaries<1) = [];
        boundaries = [0 boundaries pointrange(end)-pointrange(1)];
        for index=1:length(boundaries)-1
            tmprange = [boundaries(index)+1:boundaries(index+1)];
            EEG.data(:,tmprange) = rmbase( EEG.data(:,tmprange), length(tmprange), ...
                                               [1:length(tmprange)]);
        end;
    else
        EEG.data = rmbase( EEG.data, EEG.pnts, pointrange );    
    end;		
else
    for indc = 1:EEG.nbchan
        tmpmean  = mean(double(EEG.data(indc,pointrange,:)),2);
        EEG.data(indc,:,:) = EEG.data(indc,:,:) - repmat(tmpmean, [1 EEG.pnts 1]);
    end;
%    EEG.data = rmbase( reshape(EEG.data, EEG.nbchan, EEG.pnts*EEG.trials), EEG.pnts, pointrange );
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
