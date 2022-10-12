% POP_INTERP - interpolate data channels
%
% Usage: EEGOUT = pop_interp(EEG, badchans, method, t_range);
%
% Inputs: 
%     EEG      - EEGLAB dataset
%     badchans - [integer array] indices of channels to interpolate.
%                For instance, these channels might be bad.
%                [chanlocs structure] channel location structure containing
%                either locations of channels to interpolate or a full
%                channel structure (missing channels in the current 
%                dataset are interpolated).
%     method   - [string] method used for interpolation (default is 'spherical').
%                'invdist'/'v4' uses inverse distance on the scalp
%                'spherical' uses superfast spherical interpolation. 
%                'spacetime' uses griddata3 to interpolate both in space 
%                and time (very slow and cannot be interrupted).
%     t_range  - [integer array with just two elements] time interval of the
%                badchans which should be interpolated. First element is
%                the start time and the second element is the end time.
% Output: 
%     EEGOUT   - data set with bad electrode data replaced by
%                interpolated data
%
% Author: Arnaud Delorme, CERCO, CNRS, 2009-

% Copyright (C) Arnaud Delorme, CERCO, 2009, arno@salk.edu
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

function [EEG com] = pop_interp(EEG, bad_elec, method)

    com = '';
    if nargin < 1
        help pop_interp;
        return;
    end
    if nargin > 1 && nargin < 3
        method = 'spherical';
    end
    
    if nargin < 2
        disp('Warning: interpolation can be done on the fly in studies'); 
        disp('         this function will actually create channels in the dataset'); 
        disp('Warning: do not interpolate channels before running ICA'); 
        disp('You may define channel location to interpolate in the channel'); 
        disp('editor and declare such channels as non-data channels'); 
         
        enablenondat = 'off';
        if isfield(EEG.chaninfo, 'removedchans')
            if ~isempty(EEG.chaninfo.removedchans)
                enablenondat = 'on';
            end
        end
        if isempty(EEG.epoch)
            uilist = { { 'Style' 'text' 'string' 'What channel(s) do you want to interpolate' 'fontweight' 'bold' } ...
                       { 'style' 'text' 'string' 'none selected' 'tag' 'chanlist' } ...
                       { 'style' 'pushbutton' 'string' 'Select from removed channels' 'callback' 'pop_interp(''nondatchan'',gcbf);' 'enable' enablenondat } ...                   
                       { 'style' 'pushbutton' 'string' 'Select from data channels'    'callback' 'pop_interp(''datchan'',gcbf);' } ...                   
                       { 'style' 'pushbutton' 'string' 'Use specific channels of other dataset' 'callback' 'pop_interp(''selectchan'',gcbf);'} ...
                       { 'style' 'pushbutton' 'string' 'Use all channels from other dataset' 'callback' 'pop_interp(''uselist'',gcbf);'} ...
                       { } ...
                       { 'style' 'text'  'string' 'Interpolation method'} ...
                       { 'style' 'popupmenu'  'string' 'Spherical|Planar (slow)'  'tag' 'method' } ...
                       { } ...
                       { 'Style', 'text', 'string', 'Time range [min max] (s)', 'fontangle', fastif(length(EEG)>1, 'italic', 'normal') } ...
                       { 'Style', 'edit', 'string', '', 'enable', fastif(length(EEG)>1, 'off', 'on') } ...
                       {} { 'Style' 'text' 'string' 'Note: for group level analysis, interpolate in STUDY' } ...
                       };

            geom     = { 1 1 1 1 1 1 1 [1.1 1] 1 [1.1 1] 1   1 };
            geomvert = [ 1 1 1 1 1 1 1 1       0.5 1      0.5  1 ];
        else
            uilist = { { 'Style' 'text' 'string' 'What channel(s) do you want to interpolate' 'fontweight' 'bold' } ...
                   { 'style' 'text' 'string' 'none selected' 'tag' 'chanlist' } ...
                   { 'style' 'pushbutton' 'string' 'Select from removed channels' 'callback' 'pop_interp(''removedchans'',gcbf);' 'enable' enablenondat } ...                   
                   { 'style' 'pushbutton' 'string' 'Select from data channels'    'callback' 'pop_interp(''datchan'',gcbf);' } ...                   
                   { 'style' 'pushbutton' 'string' 'Use specific channels of other dataset' 'callback' 'pop_interp(''selectchan'',gcbf);'} ...
                   { 'style' 'pushbutton' 'string' 'Use all channels from other dataset' 'callback' 'pop_interp(''uselist'',gcbf);'} ...
                   { } ...
                   { 'style' 'text'  'string' 'Interpolation method'} ...
                   { 'style' 'popupmenu'  'string' 'Spherical|Planar (slow)'  'tag' 'method' } ...
                   {} { 'Style' 'text' 'string' 'Note: for group level analysis, interpolate in STUDY' } ...
                   };
               
            geom     = { 1 1 1 1 1 1 1 [1.1 1] 1   1 };
            geomvert = [ 1 1 1 1 1 1 1 1       0.5 1 ];
        end
        [res, userdata, ~, restag ] = inputgui( 'uilist', uilist, 'title', 'Interpolate channel(s) -- pop_interp()', 'geometry', geom, 'geomvert', geomvert, 'helpcom', 'pophelp(''pop_interp'')');
        if isempty(res) || isempty(userdata), return; end
        
        if restag.method == 1
             method = 'spherical';
        else method = 'invdist';
        end
        bad_elec = userdata.chans;
        if numel(res) > 1
            t_range = res{2};
        else
            t_range = '';
        end
        if isempty(t_range)
            t_range = [EEG.xmin EEG.xmax];
        else
            t_range = eval( [ '[' t_range ']' ] );
        end
        if size(t_range,2) ~= 2
            error('Time/point range must contain 2 columns exactly');
        end
        if floor(max(t_range)) > EEG.xmax 
            error('Time/point range exceed upper data limits');
        end
        if min(t_range) < EEG.xmin
            error('Time/point range exceed lower data limits');
        end
        
        com = sprintf('EEG = pop_interp(EEG, %s, ''%s'');', userdata.chanstr, method);
        if ~isempty(findstr('removedchans', userdata.chanstr))
            eval( [ userdata.chanstr '=[];' ] );
        end
        
    elseif ischar(EEG)
        command = EEG;
        clear EEG;
        fig = bad_elec;
        userdata = get(fig, 'userdata');
        
        if strcmpi(command, 'removedchans')
            global EEG;
            tmpchaninfo = EEG.chaninfo;
            [chanlisttmp chanliststr] = pop_chansel( { tmpchaninfo.removedchans.labels } );
            if ~isempty(chanlisttmp),
                userdata.chans   = EEG.chaninfo.removedchans(chanlisttmp);
                userdata.chanstr = [ 'EEG.chaninfo.removedchans([' num2str(chanlisttmp) '])' ];
                set(fig, 'userdata', userdata);
                set(findobj(fig, 'tag', 'chanlist'), 'string', chanliststr);
            end
        elseif strcmpi(command, 'datchan')
            global EEG;
            tmpchaninfo = EEG.chanlocs;
            [chanlisttmp chanliststr] = pop_chansel( { tmpchaninfo.labels } );
            if ~isempty(chanlisttmp),
                userdata.chans   = chanlisttmp;
                userdata.chanstr = [ '[' num2str(chanlisttmp) ']' ];
                set(fig, 'userdata', userdata);
                set(findobj(fig, 'tag', 'chanlist'), 'string', chanliststr);
            end
        else
            global ALLEEG EEG;
            tmpanswer = inputdlg2({ 'Dataset index' }, 'Choose dataset', 1, { '' });
            if ~isempty(tmpanswer),
                tmpanswernum = round(str2num(tmpanswer{1}));
                if ~isempty(tmpanswernum),
                    if tmpanswernum > 0 && tmpanswernum <= length(ALLEEG),
                        TMPEEG = ALLEEG(tmpanswernum);
                        
                        tmpchans1 = TMPEEG.chanlocs;
                        if strcmpi(command, 'selectchan')
                            chanlist = pop_chansel( { tmpchans1.labels } );
                        else
                            chanlist = 1:length(TMPEEG.chanlocs); % use all channels
                        end
                        
                        % look at what new channels are selected
                        tmpchans2 = EEG.chanlocs;
                        [tmpchanlist chaninds] = setdiff_bc( { tmpchans1(chanlist).labels }, { tmpchans2.labels } );
                        if ~isempty(tmpchanlist),
                            if length(chanlist) == length(TMPEEG.chanlocs)
                                userdata.chans   = TMPEEG.chanlocs;
                                userdata.chanstr = [ 'ALLEEG(' tmpanswer{1} ').chanlocs' ];
                            else
                                userdata.chans   = TMPEEG.chanlocs(chanlist(sort(chaninds)));
                                userdata.chanstr = [ 'ALLEEG(' tmpanswer{1} ').chanlocs([' num2str(chanlist(sort(chaninds))) '])' ];
                            end
                            set(fig, 'userdata', userdata);
                            tmpchanlist(2,:) = { ' ' };
                            set(findobj(gcbf, 'tag', 'chanlist'), 'string', [ tmpchanlist{:} ]);
                        else
                            warndlg2('No new channels selected');
                        end
                    else
                        warndlg2('Wrong index');
                    end
                end
            end
        end
        return;
    end
    
    % remove from removedchans if interpolated
    if isfield(EEG.chaninfo, 'removedchans') && ~isempty(EEG.chaninfo.removedchans) && isstruct(bad_elec)
        for iChan = 1:length(bad_elec)
            ind = strmatch(lower(bad_elec(iChan).labels), lower({EEG.chaninfo.removedchans.labels}), 'exact');
            if ~isempty(ind)
                EEG.chaninfo.removedchans(ind) = [];
            end
        end
    end

    EEG = eeg_interp(EEG, bad_elec, method, t_range);
    
    
