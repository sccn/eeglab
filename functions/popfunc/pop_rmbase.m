% pop_rmbase() - remove channel baseline means from an epoched or 
%                continuous EEG dataset. Calls rmbase().
% Usage:
%   >> OUTEEG = pop_rmbase( EEG ); % pop up an interactive arg entry window
%   >> OUTEEG = pop_rmbase( EEG, timerange, pointrange, chanlist); % call rmbase()
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
%   chanlist   - [cell]          List of channels. Default is all.
%
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

function [EEG, com] = pop_rmbase( EEG, timerange, pointrange, chanlist)

com ='';
if nargin < 1
    help pop_rmbase;
    return;
end
if isempty(EEG(1).data)
    disp('pop_rmbase(): cannot remove baseline of an empty dataset'); return;
end    
if nargin < 1
	help pop_rmbase;
	return;
end
if nargin < 4 || isempty(chanlist)
    chanlist = 1:EEG(1).nbchan;
end
if nargin < 2, timerange = [];  end
if nargin < 3, pointrange = []; end
if nargin < 2 
    % popup window parameters
    % -----------------------
    defaultbase = [num2str(EEG(1).times(1)) ' 0'];
    if EEG(1).times(1) >= 0
        defaultbase = '[ ]';
    end
    if EEG(1).trials > 1
        uilist = { { 'style' 'text' 'string' 'Baseline latency range ([min max] in ms) ([ ] = whole epoch):' } ...
            { 'style' 'edit' 'string'  defaultbase } ...
            { 'style' 'text' 'string' 'Or remove baseline points vector (ex:1:56):' } ...
            { 'style' 'edit' 'string' '' } ...
            { 'style' 'text' 'string' 'Note: press Cancel if you do not want to remove the baseline' } };
        uigeom = { [3 1] [3 1] [1] };
    else
        uilist = { { 'style' 'text' 'string' 'Removing the mean of each data channel (press cancel to skip)' } };
        uigeom = { [1] };
    end
    
    % add channel selection
    % --------------
    if ~isempty(EEG(1).chanlocs)
        tmpchanlocs = EEG(1).chanlocs;        
    else
        tmpchanlocs = [];
        for index = 1:EEG(1).nbchan
            tmpchanlocs(index).labels = int2str(index);
            tmpchanlocs(index).type = '';
        end
    end
    cb_type = 'pop_chansel(get(gcbf, ''userdata''), ''field'', ''type'',   ''handle'', findobj(''parent'', gcbf, ''tag'', ''chantypes''));';
    cb_chan = 'pop_chansel(get(gcbf, ''userdata''), ''field'', ''labels'', ''handle'', findobj(''parent'', gcbf, ''tag'', ''channels''));';
    enableChan = fastif(length(EEG) > 1, 'off', 'on') ;
    uilist = { uilist{:} ...
        { 'style' 'text'       'string' 'Channel type(s)' 'enable' enableChan} ...
        { 'style' 'edit'       'string' '' 'tag' 'chantypes' 'enable' enableChan}  ...
        { 'style' 'pushbutton' 'string' '...'  'callback' cb_type  'enable' enableChan} ...
        { 'style' 'text'       'string' 'OR channel(s) (default all)'  'enable' enableChan} ...
        { 'style' 'edit'       'string' '' 'tag' 'channels'  'enable' enableChan}  ...
        { 'style' 'pushbutton' 'string' '...' 'callback' cb_chan  'enable' enableChan}
        };
    uigeom = { uigeom{:} [2 1.5 0.5] [2 1.5 0.5] };
    [result, usrdat, sres2, sres] = inputgui( 'uilist', uilist, 'geometry', uigeom, 'title', 'Baseline removal - pop_rmbase()', 'helpcom', 'pophelp(''pop_rmbase'');', 'userdata', tmpchanlocs);
    if isempty(result), return; end
    
    % decode parameters
    % -----------------
    if EEG(1).trials > 1
        if numel(result) < 2 || ((isempty(result{1}) || strcmp(result{1},'[]') ) ...
                && (isempty(result{2}) || strcmp(result{2},'[]')))
            timerange = [ EEG(1).times(1) EEG(1).times(end) ]; % whole epoch latency range
            pointrange = [];
            fprintf('pop_rmbase(): using whole epoch as baseline.\n');
        else
            timerange  = eval( [ '[' result{1} ']' ] );
            pointrange = eval( [ '[' result{2} ']' ] );
        end
    end
    if ~isempty(sres.chantypes)
        chanlist = eeg_decodechan(EEG.chanlocs, parsetxt(sres.chantype), 'type');
    elseif ~isempty(sres.channels)
        chanlist = eeg_decodechan(EEG(1).chanlocs, sres.channels);
    end
    
end

% process multiple datasets
% -------------------------
if length(EEG) > 1
    if nargin < 2
        [ EEG, com ] = eeg_eval( 'pop_rmbase', EEG, 'warning', 'on', 'params', { timerange pointrange } );
    else
        [ EEG, com ] = eeg_eval( 'pop_rmbase', EEG, 'params', { timerange pointrange } );
    end
    return;
end

flag_timerange = 1; % provide timerange as input (so use in history)
if ~isempty(timerange)
    if timerange(1) < EEG.times(1) || timerange(end) > EEG.times(end)
        error('pop_rmbase(): Bad time range');
    end
    pointrange = find( EEG.times >= timerange(1) & EEG.times <= timerange(2));
elseif ~isempty(pointrange)
    flag_timerange = 0;
    if pointrange(1) < 1, pointrange(1) = 1; end
    if pointrange(end) > EEG.pnts, pointrange(end) = EEG.pnts; end
else
    pointrange = [1:EEG.pnts];
end

%
% Respect excised data boundaries if continuous data
% ---------------------------------------------------
fprintf('pop_rmbase(): Removing baseline...\n');
if EEG.trials == 1 && ~isempty(EEG.event) ...
                     && isfield(EEG.event, 'type') ...
                        && ischar(EEG.event(1).type)
    tmpevent = EEG.event;
	boundaries = strmatch('boundary', {tmpevent.type});
	if ~isempty(boundaries) % this is crashing
        fprintf('Pop_rmbase(): finding continuous data discontinuities\n');
        boundaries = round([ tmpevent(boundaries).latency ] -0.5-pointrange(1)+1);
        boundaries(boundaries>=pointrange(end)-pointrange(1)) = [];
        boundaries(boundaries<1) = [];
        boundaries = [0 boundaries pointrange(end)-pointrange(1)+1];
        for index=1:length(boundaries)-1
            tmprange = [boundaries(index)+1:boundaries(index+1)];
            if length(tmprange) > 1
                EEG.data(chanlist,tmprange) = rmbase( EEG.data(:,tmprange), length(tmprange), ...
                                                   [1:length(tmprange)]);
            elseif length(tmprange) == 1
                EEG.data(chanlist,tmprange) = 0;
            end
        end
    else
        EEG.data(chanlist,:) = rmbase( EEG.data(chanlist,:), EEG.pnts, pointrange );    
    end
else
    for indc = chanlist
        tmpmean  = mean(double(EEG.data(indc,pointrange,:)),2);
        EEG.data(indc,:,:) = EEG.data(indc,:,:) - repmat(tmpmean, [1 EEG.pnts 1]);
    end
%    EEG.data = rmbase( reshape(EEG.data, EEG.nbchan, EEG.pnts*EEG.trials), EEG.pnts, pointrange );
end

EEG.data = reshape( EEG.data, EEG.nbchan, EEG.pnts, EEG.trials);
EEG.icaact = [];

if isequal(chanlist, [1:EEG.nbchan]), chanlist = []; end
if flag_timerange, pointrange = []; else, timerange = []; end
if isempty(chanlist)
    com = sprintf('EEG = pop_rmbase( EEG, %s);',vararg2str({timerange pointrange}));
else
    com = sprintf('EEG = pop_rmbase( EEG, %s);',vararg2str({timerange pointrange chanlist}));
end
return;
