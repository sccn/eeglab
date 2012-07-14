% pop_interp() - interpolate data channels
%
% Usage: EEGOUT = pop_interp(EEG, badchans, method);
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
%                and time (very slow and cannot be interupted).
% Output: 
%     EEGOUT   - data set with bad electrode data replaced by
%                interpolated data
%
% Author: Arnaud Delorme, CERCO, CNRS, 2009-

% Copyright (C) Arnaud Delorme, CERCO, 2009, arno@salk.edu
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

function [EEG com] = pop_interp(EEG, bad_elec, method)

    com = '';
    if nargin < 1
        help pop_interp;
        return;
    end;
    
    if nargin < 2
        disp('Warning: interpolation can be done on the fly in studies'); 
        disp('         this function will actually create channels in the dataset'); 
        disp('Warning: do not interpolate channels before running ICA'); 
        disp('You may define channel location to interpolate in the channel'); 
        disp('editor and declare such channels as non-data channels'); 
         
        enablenondat = 'off';
        if isfield(EEG.chaninfo, 'nodatchans')
            if ~isempty(EEG.chaninfo.nodatchans)
                enablenondat = 'on';
            end;
        end;
        cb_nondat   = [ 'tmpchaninfo = EEG.chaninfo; [chanlisttmp chanliststr] = pop_chansel( { tmpchaninfo.nodatchans.labels } );' ...
                      'if ~isempty(chanlisttmp),' ...
                      '   set(gcbf, ''userdata'', EEG.chaninfo.nodatchans(chanlisttmp));' ...
                      '   set(findobj(gcbf, ''tag'', ''chanlist''), ''string'', chanliststr);' ...                      
                      'end;' ...
                      'clear chanlisttmp chanliststr tmpchaninfo;' ];
                      
        cb_otherdat = [ 'tmpanswer = inputdlg2({ ''Dataset index'' }, ''Choose dataset'', 1, { '''' });' ...
                      'if ~isempty(tmpanswer),' ...
                      '   tmpanswernum = round(str2num(tmpanswer{1}));' ...
                      '   if ~isempty(tmpanswernum),' ...
                      '       if tmpanswernum > 0 & tmpanswernum < length(ALLEEG),' ...
                      '           tmpchans1 = ALLEEG(tmpanswernum).chanlocs;' ...
                      '           tmpchans2 = EEG.chanlocs;' ...
                      '           tmpchanlist = setdiff( { tmpchans1.labels }, { tmpchans2.labels } );' ...
                      '           if ~isempty(tmpchanlist),' ...
                      '               set(gcbf, ''userdata'', ALLEEG(tmpanswernum).chanlocs);' ...
                      '               tmpchanlist(2,:) = { '' '' };' ...
                      '               set(findobj(gcbf, ''tag'', ''chanlist''), ''string'', [ tmpchanlist{:} ]);' ...
                      '           else,' ...
                      '               warndlg2(''No new channels in this dataset'');' ...
                      '           end;' ...
                      '       else,' ...
                      '           warndlg2(''Wrong index'');' ...
                      '       end;' ...
                      '   end;' ...
                      'end;' ...
                      'clear tmpanswer tmpanswernum;' ];
                  
        uilist = { { 'Style' 'text' 'string' 'What channel(s) to interpolate' 'fontweight' 'bold' } ...
                   { 'style' 'text' 'string' 'none' 'tag' 'chanlist' } ...
                   { 'style' 'pushbutton' 'string' 'Select from non-data channels' 'callback' 'pop_interp(''nondatchan'',gcbf);' 'enable' enablenondat } ...                   
                   { 'style' 'pushbutton' 'string' 'Select from other dataset' 'callback' 'pop_interp(''selectchan'',gcbf);'} ...
                   { 'style' 'pushbutton' 'string' 'Use list of other dataset' 'callback' 'pop_interp(''uselist'',gcbf);'} ...
                   { } ...
                   { 'style' 'text'  'string' 'Interpolation method'} ...
                   { 'style' 'popupmenu'  'string' 'Spherical|Planar (slow)'  'tag' 'method' } ...
                   };
               
        geom = { 1 1 1 1 1 1 [1.8 1] };
        [res userdata tmp restag ] = inputgui( 'uilist', uilist, 'title', 'Interpolate channel(s) -- pop_interp()', 'geometry', geom, 'helpcom', 'pophelp(''pop_interp'')');
        if isempty(res) | isempty(userdata), return; end;
        
        if restag.method == 1
             method = 'spherical';
        else method = 'invdist';
        end;
        bad_elec = userdata.chans;
        
        com = sprintf('EEG = pop_interp(EEG, %s, ''%s'');', userdata.chanstr, method);
        if ~isempty(findstr('nodatchans', userdata.chanstr))
            eval( [ userdata.chanstr '=[];' ] );
        end;
        
    elseif isstr(EEG)
        command = EEG;
        clear EEG;
        fig = bad_elec;
        userdata = get(fig, 'userdata');
        
        if strcmpi(command, 'nondatchan')
            global EEG;
            tmpchaninfo = EEG.chaninfo;
            [chanlisttmp chanliststr] = pop_chansel( { tmpchaninfo.nodatchans.labels } );
            if ~isempty(chanlisttmp),
                userdata.chans   = EEG.chaninfo.nodatchans();
                userdata.chanstr = [ 'EEG.chaninfo.nodatchans([' num2str(chanlisttmp) '])' ];
                set(fig, 'userdata', userdata);
                set(findobj(fig, 'tag', 'chanlist'), 'string', chanliststr);
            end;
        else
            global ALLEEG EEG;
            tmpanswer = inputdlg2({ 'Dataset index' }, 'Choose dataset', 1, { '' });
            if ~isempty(tmpanswer),
                tmpanswernum = round(str2num(tmpanswer{1}));
                if ~isempty(tmpanswernum),
                    if tmpanswernum > 0 & tmpanswernum < length(ALLEEG),
                        TMPEEG = ALLEEG(tmpanswernum);
                        
                        tmpchans1 = TMPEEG.chanlocs;
                        if strcmpi(command, 'selectchan')
                            chanlist = pop_chansel( { tmpchans1.labels } );
                        else
                            chanlist = 1:length(TMPEEG.chanlocs); % use all channels
                        end
                        
                        % look at what new channels are selected
                        tmpchans2 = EEG.chanlocs;
                        [tmpchanlist chaninds] = setdiff( { tmpchans1(chanlist).labels }, { tmpchans2.labels } );
                        if ~isempty(tmpchanlist),
                            if length(chanlist) == length(TMPEEG.chanlocs)
                                userdata.chans   = TMPEEG.chanlocs;
                                userdata.chanstr = [ 'ALLEEG(' tmpanswer{1} ').chanlocs' ];
                            else
                                userdata.chans   = TMPEEG.chanlocs(chanlist(sort(chaninds)));
                                userdata.chanstr = [ 'ALLEEG(' tmpanswer{1} ').chanlocs([' num2str(chanlist(sort(chaninds))) '])' ];
                            end;
                            set(fig, 'userdata', userdata);
                            tmpchanlist(2,:) = { ' ' };
                            set(findobj(gcbf, 'tag', 'chanlist'), 'string', [ tmpchanlist{:} ]);
                        else
                            warndlg2('No new channels selected');
                        end;
                    else
                        warndlg2('Wrong index');
                    end;
                end;
            end;
        end;
        return;
    end;
    
    EEG = eeg_interp(EEG, bad_elec, method);
    
    
