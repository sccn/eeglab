% pop_reref() - Convert an EEG dataset to average reference or to a
%               new common reference channel (or channels). Calls reref().
% Usage:
%       >> EEGOUT = pop_reref( EEG ); % pop up interactive window
%       >> EEGOUT = pop_reref( EEG, ref, 'key', 'val' ...);
%
% Graphic interface:
%   "Compute average reference" - [edit box] Checking this box (for 'yes') is 
%                 the same as giving an empty value for the commandline 'ref' 
%                 argument. Unchecked, the data are transformed to common reference.
%   "Re-reference data to channel(s)" - [checkbox] Checking this option
%                 automatically unchecks the checkbox above, allowing reference 
%                 channel indices to be entered in the text edit box to its right
%                 (No commandline equivalent).
%   "Interpolate removed channel(s)"  - [checkbox] Enable the interpolation
%                 of removed channels to compute the re-referencing.
%                 Checking this option is equivalent to use the command
%                 line option 'interpchan' with argument '[]'
%                 (see optional input 'interpchan')
%   "Retain old reference channels in data" - [checkbox] When re-referencing the 
%                 data, checking this checkbox includes the data for the 
%                 previous reference channel.
%   "Exclude channel indices (EMG, EOG)" - [edit box] exclude the given
%                 channel indices from rereferencing.
%   "Add current reference channel back to the data" - [edit box] When 
%                 re-referencing the data, checking this checkbox
%                 reconstitutes the data for the previous reference
%                 channel. If the location for this channel or the channel is
%                 not present, it first needs to be defined in the channel
%                 editor (additional channel defined as reference and not 
%                 associated with data).
% Inputs:
%   EEG         - input dataset
%   ref         - reference: []             = convert to average reference
%                            [int vector]   = new reference electrode number(s)
%                            'Cz'           = string
%                            { 'P09' 'P10 } = cell array of strings
% Optional inputs:
%   'interpchan'  - [channel location structure | integer array | [] | 'off']
%                   Channels to interpolate prior to re-referencing the data. If [],
%                   channels will be found by comparing all the channels (type = EEG) 
%                   in the current EEG.chanlocs structure against EEG.urchanlocs. A channel 
%                   location  structure of the channels to be interpolated can be provided 
%                   as an input, as well as the index of the channels into the EEG.urchanlocs.
%                   Default:'off'
%   'exclude'     - [integer array] List of channels to exclude. Default: none.
%   'keepref'     - ['on'|'off'] keep the reference channel. Default: 'off'.
%   'refloc'      - [structure] Previous reference channel structure. Default: none.
%
% Outputs:
%   EEGOUT      - re-referenced output dataset
%
% Notes:
%                 For other options, call reref() directly. See >> help reref
%
% Author: Arnaud Delorme, CNL / Salk Institute, 12 Nov 2002
%
% See also: reref(), eeglab()

% Copyright (C) 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [EEG, com] = pop_reref( EEG, ref, varargin);

com = '';
if nargin < 1
   help pop_reref;
   return;
end;   
if isempty(EEG(1).data)
    error('Pop_reref: cannot process empty data');
end

% gui inputs
% ----------
if nargin < 2
    
    % find initial reference
    % ----------------------
    if length(EEG(1).chanlocs) == EEG(1).nbchan+1
        includeref = 1;
    end
    
    geometry = { [1] [1] [1.8 1 0.3] [1] [1] [1] [1.8 1 0.3] [1.8 1 0.3] };
    cb_setref = [ 'set(findobj(''parent'', gcbf, ''tag'', ''refbr'')    , ''enable'', ''on'');' ...
                  'set(findobj(''parent'', gcbf, ''tag'', ''reref'')    , ''enable'', ''on'');' ...
                  'set(findobj(''parent'', gcbf, ''tag'', ''keepref'')  , ''enable'', ''on'');' ];
    cb_setave = [ 'set(findobj(''parent'', gcbf, ''tag'', ''refbr'')    , ''enable'', ''off'');' ...
                  'set(findobj(''parent'', gcbf, ''tag'', ''reref'')    , ''enable'', ''off'');' ...
                  'set(findobj(''parent'', gcbf, ''tag'', ''keepref'')  , ''enable'', ''off'', ''value'', 0);' ];
    cb_averef = [ 'set(findobj(''parent'', gcbf, ''tag'', ''rerefstr'') , ''value'', ~get(gcbo, ''value''));' ...
                  'if get(gcbo, ''value''),' cb_setave ...
                  'else,'                    cb_setref ...
                  'end;' ];
    cb_ref    = [ 'set(findobj(''parent'', gcbf, ''tag'', ''ave'')      , ''value'', ~get(gcbo, ''value''));' ...
                  'if get(gcbo, ''value''),' cb_setref ...
                  'else,'                    cb_setave ...
                  'end;' ];
    cb_chansel1 = 'tmpEEG = get(gcbf, ''userdata''); tmpchanlocs = tmpEEG(1).chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels}, ''withindex'', ''on''); set(findobj(gcbf, ''tag'', ''reref''   ), ''string'',tmpval); clear tmpEEG tmpchanlocs tmp tmpval';
    cb_chansel2 = 'tmpEEG = get(gcbf, ''userdata''); tmpchanlocs = tmpEEG(1).chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels}, ''withindex'', ''on''); set(findobj(gcbf, ''tag'', ''exclude'' ), ''string'',tmpval); clear tmpEEG tmpchanlocs tmp tmpval';
    cb_chansel3 = [ 'tmpEEG = get(gcbf, ''userdata''); if ~isfield(tmpEEG(1).chaninfo, ''nodatchans''), ' ...
                    '   warndlg2(''There are no Reference channel defined, add it using the channel location editor'');' ...
                    'elseif isempty(tmpEEG(1).chaninfo.nodatchans),' ...
                    '   warndlg2(''There are no Reference channel defined, add it using the channel location editor'');' ...
                    'elseif isfield(tmpEEG(1).chaninfo.nodatchans, ''type''),' ...
                    '   fidType = ismember(cellfun(@char, {  tmpEEG(1).chaninfo.nodatchans.type}, ''UniformOutput'', false), ''FID'');' ...
                    '   if sum(fidType == 0) == 0,' ...
                    '      warndlg2(''There are no Reference channel defined, add it using the channel location editor'');' ...
                    '   else,' ...
                    '      tmpchaninfo = tmpEEG(1).chaninfo; [tmp tmpval] = pop_chansel({tmpchaninfo.nodatchans(~fidType).labels}, ''withindex'', ''on'');' ...
                    '      set(findobj(gcbf, ''tag'', ''refloc''  ), ''string'',tmpval);' ...
                    '   end;' ...
                    'end;' ...
                    'clear tmpEEG tmpchanlocs tmp tmpval;' ];
    if isempty(EEG(1).chanlocs), cb_chansel1 = ''; cb_chansel2 = ''; cb_chansel3 = ''; end
    
    % find current reference (= reference most used)
    % ----------------------------------------------
    if isfield(EEG(1).chanlocs, 'ref')
        tmpchanlocs = EEG(1).chanlocs;
        [curref,~,allinds] = unique_bc( { tmpchanlocs.ref });
        maxind = 1;
        for ind = unique_bc(allinds)
            if length(find(allinds == ind)) > length(find(allinds == maxind))
                maxind = ind;
            end
        end
        curref = curref{maxind};
        if isempty(curref), curref = 'unknown'; end
    else
        curref = 'unknown';
    end
    
    uilist = { { 'style' 'text' 'string' [ 'Current data reference state is: ' curref] } ...
               ...
               { 'style' 'checkbox' 'tag' 'ave'   'value' 1 'string' 'Compute average reference' 'callback' cb_averef } ...
               ...
               { 'style' 'checkbox' 'tag' 'rerefstr' 'value' 0 'string' 'Re-reference data to channel(s):' 'callback'  cb_ref } ...
               { 'style' 'edit' 'tag' 'reref' 'string' '' 'enable' 'off' } ...
               { 'style' 'pushbutton' 'string' '...' 'callback' cb_chansel1 'enable' 'off' 'tag' 'refbr' } ...
               { 'style' 'checkbox' 'tag' 'interp' 'value' 0 'string' 'Interpolate removed channel(s)'} ...
               ...
               {} ...
               ...
               { 'style' 'checkbox' 'value' 0 'enable' 'off' 'tag' 'keepref' 'string' 'Retain ref. channel(s) in data (will be flat for single-channel ref.)' } ...
               ...
               { 'style' 'text' 'string' 'Exclude channel indices (EMG, EOG)' } ...
               { 'style' 'edit' 'tag' 'exclude' 'string' '' } ...
               { 'style' 'pushbutton' 'string' '...' 'callback' cb_chansel2 } ...
               ...
               { 'style' 'text' 'tag' 'reflocstr' 'string' 'Add old ref. channel back to the data' } ...
               { 'style' 'edit' 'tag' 'refloc' 'string' '' } ...
               { 'style' 'pushbutton' 'string' '...' 'callback' cb_chansel3 } };
    
    [result,~,~,restag] = inputgui('geometry', geometry, 'uilist', uilist, 'helpcom', 'pophelp(''pop_reref'')', 'title', 'pop_reref - average reference or re-reference data', 'userdata', EEG);
    if isempty(result), return; end

    % decode inputs
    % -------------
    options = {};
    if ~isempty(restag.refloc)
         try
             tmpchaninfo = EEG(1).chaninfo;
             tmpallchans = lower({ tmpchaninfo.nodatchans.labels });
             allelecs = parsetxt(lower(restag.refloc));
             chanind  = [];
             for iElec = 1:length(allelecs)
                 chanind = [chanind strmatch( allelecs{iElec}, tmpallchans, 'exact') ];
             end
             options = { options{:} 'refloc' EEG(1).chaninfo.nodatchans(chanind) }; 
         catch, disp('Error with old reference: ignoring it');
         end
    end
    if ~isempty(restag.exclude), options = { options{:} 'exclude' eeg_chaninds(EEG, restag.exclude) }; end
    if restag.keepref,           options = { options{:} 'keepref' 'on' }; end
    if restag.ave,               ref = []; end
    if restag.rerefstr           
        if isempty(restag.reref)
            warndlg2('Aborting: you must enter one or more reference channels'); 
            return;
        else
            ref = eeg_chaninds(EEG, restag.reref); 
        end
    end
    if restag.interp == 1, options = { options{:} 'interpchan' [] }; end
else
    options = varargin;
end
if ischar(ref), ref = { ref }; end

% process multiple datasets
% -------------------------
if length(EEG) > 1
    if nargin < 2
        [ EEG, com ] = eeg_eval( 'pop_reref', EEG, 'warning', 'on', 'params', {ref options{:} } );
    else
        [ EEG, com ] = eeg_eval( 'pop_reref', EEG, 'params', {ref options{:} } );
    end
    return;
end

orichanlocs = EEG.chanlocs;
orinbchan   = EEG.nbchan;
if iscell(ref), ref = eeg_chaninds(EEG, ref); end
optionscall = options;

g = struct(optionscall{:});
if ~isfield(g, 'exclude'),       g.exclude       = [];    end
if ~isfield(g, 'keepref'),       g.keepref       = 'off'; end
if ~isfield(g, 'refloc') ,       g.refloc        = [];    end
if ~isfield(g, 'interpchan') ,   g.interpchan    = 'off'; end
if ~isfield(g, 'addrefchannel'), g.addrefchannel = 0;     end
if ~isfield(g, 'enforcetype'),   g.enforcetype   = 0;     end
%--- Interpolation code START
interpflag = 0;
if ~isequal('off', g.interpchan )
    
    % Case no channel provided, inferring them from urchanlocs field
    if isempty(g.interpchan) 
        if isfield(EEG.chaninfo, 'removedchans')
            chanlocs2interp = EEG.chaninfo.removedchans;
            if ~isempty(chanlocs2interp)
                interpflag = 1;
            end
        else
            try
                urchantype  = {EEG.urchanlocs.type};
                chanloctype = {EEG.chanlocs.type};
                if  any(cellfun(@isempty,urchantype)) || any(cellfun(@isempty,chanloctype))
                    eegtypeindx0 = [1:length(EEG.urchanlocs)]';
                    eegtypeindx1 = [1:length(EEG.chanlocs)]';
                    % Excluding fiducials if exist
                    try
                        indxfid_urch = find(strcmpi({'fid'},urchantype));
                        indxfid_ch   = find(strcmpi({'fid'},chanloctype));
                        if ~isempty(indxfid_urch), eegtypeindx0(indxfid_urch) = []; end
                        if ~isempty(indxfid_ch),   eegtypeindx1(indxfid_ch)   = []; end
                    catch
                        fprintf('pop_reref message: Unable to find fiducials...\n');
                    end
                else   
                    eegtypeindx0 = strmatch('EEG',urchantype);
                    eegtypeindx1 = strmatch('EEG',chanloctype);    
                end     
            catch
                fprintf(2,'pop_reref error: Unable to check for deleted channels. Missing field ''type'' in channel location \n');
                return;
            end

            if isempty(eegtypeindx0) || isempty(eegtypeindx1)
                fprintf(2,'pop_reref error: Unable to get channel type from this data. Check field ''type'' on EEG.urchanlocs or EEG.chanlocs. \n');
                return;
            end

            chan2interp = setdiff_bc({EEG.urchanlocs(eegtypeindx0).labels}, {EEG.chanlocs(eegtypeindx1).labels});       
            if isempty(chan2interp)
                fprintf('pop_reref message: No removed channel found. Halting interpolation and moving forward...\n');
            else
               chan2interpindx = find(cell2mat(cellfun(@(x) ismember(x, chan2interp), {EEG.urchanlocs.labels}, 'UniformOutput', 0)));  

               % Checking validity of channels selected for interpolation by assessing X coordinate
               for ichan = 1:length(chan2interpindx)
                   validchan(ichan) = ~isempty(EEG.urchanlocs(chan2interpindx(ichan)).X);
               end
               if all(validchan == 0)
                   fprintf('pop_reref message: Invalid channel(s) for interpolation. Halting interpolation and moving forward...\n');
               else
                   chanlocs2interp =  EEG.urchanlocs(chan2interpindx(find(validchan)));
                   interpflag = 1;
               end
            end
        end
    % Case where channel loc structure is provided    
    elseif isstruct(g.interpchan) 
        chanlocs2interp = g.interpchan;
        interpflag = 1;
    
    % Case where channel index is provided    
    elseif isreal(g.interpchan)
        chanlocs2interp = EEG.urchanlocs(g.interpchan);
        interpflag = 1;
        
    % invalid case    
    else
        error('pop_reref error: Invalid arguments for option ''interpchan'' ');
    end
    
    if interpflag
        EEG = pop_interp(EEG, chanlocs2interp, 'spherical');
        interpindx = find(cell2mat(cellfun(@(x) ismember(x, {chanlocs2interp.labels}), {EEG.chanlocs.labels}, 'UniformOutput', 0))); 
        % EEG.icaweights = [];
    end
end
%--- interpolation code END

% include channel location file
% -----------------------------
if ~isempty(EEG.chanlocs)
    optionscall = { optionscall{:} 'elocs' EEG.chanlocs }; 
end;    

fprintf('Re-referencing data\n');
[EEG.data EEG.chanlocs refchan ] = reref(EEG.data, ref, optionscall{:});

% If interpolation was done... then remove channels
if interpflag
    EEG = pop_select(EEG, 'nochannel', interpindx);
end

nchans = EEG.nbchan; % retrieve number of channels for ICA bussines

% deal with reference
% -------------------
if ~isempty(refchan)
    if ~isfield(EEG.chaninfo, 'nodatchans')
        EEG.chaninfo.nodatchans = refchan;
    elseif isempty(EEG.chaninfo.nodatchans)
        EEG.chaninfo.nodatchans = refchan;
    else
        allf = fieldnames(refchan);
        n    = length(EEG.chaninfo.nodatchans);
        for iRef = 1:length(refchan)
            for ind = 1:length(allf)
                EEG.chaninfo.nodatchans = setfield(EEG.chaninfo.nodatchans, { n+iRef }, ...
                    allf{ind}, getfield(refchan(iRef), allf{ind}));
            end
        end
    end
end
if ~isempty(g.refloc) 
    if isfield(EEG.chaninfo, 'nodatchans') && ~isempty(EEG.chaninfo.nodatchans)
        allinds = [];
        tmpchaninfo = EEG.chaninfo;
        for iElec = 1:length(g.refloc)
            if isempty(tmpchaninfo) || isempty(tmpchaninfo.nodatchans)
                error('Missing reference channel information. Edit channels and add reference first.');
            end
            allinds = [allinds strmatch( g.refloc(iElec).labels, { tmpchaninfo.nodatchans.labels }) ];
        end
        EEG.chaninfo.nodatchans(allinds) = [];
    else
        error('Missing reference channel information. Edit channels and add reference first.');
    end
end
    
% legacy EEG.ref field
% --------------------
if isfield(EEG, 'ref')
    if strcmpi(EEG.ref, 'common') && isempty(ref)
        EEG.ref = 'average';
    elseif strcmpi(EEG.ref, 'average') && ~isempty(ref)
        EEG.ref = 'common';
    end
end

EEG.nbchan = size(EEG.data,1);
EEG = eeg_checkset(EEG);

% include ICA or not
% ------------------
if ~isempty(EEG.icaweights)
    
    if ~isempty(intersect(EEG.icachansind, g.exclude))
        disp('Warning: some channels used for ICA were excluded from referencing');
        disp('         the ICA decomposition has been removed');
        EEG.icaweights = [];
        EEG.icasphere  = [];
    elseif length(EEG.icachansind) ~= nchans - length(g.exclude)
        disp('Error: some channels not used for ICA decomposition are used for rereferencing');
        disp('       the ICA decomposition has been removed');
        EEG.icaweights = [];
        EEG.icasphere  = [];
    else
        fprintf('Re-referencing ICA matrix\n');
        if isempty(orichanlocs)
            error('Cannot re-reference ICA decomposition without channel locations')
        end
        
        newICAchaninds = zeros(orinbchan, size(EEG.icawinv,2));
        newICAchaninds(EEG.icachansind,:) = EEG.icawinv;
        
        [newICAchaninds, newchanlocs] = reref(newICAchaninds, ref, optionscall{:});
        
        % convert channel indices in icachanlocs (uses channel labels)
        % ------------------------------------------------------------
        icachansind = EEG.icachansind;
        rminds      = [1:size(newICAchaninds,1)];
        for i=length(icachansind):-1:1
            oldLabel    = orichanlocs(icachansind(i)).labels;
            newLabelPos = strmatch(oldLabel, { newchanlocs.labels }, 'exact');
            
            if length(newLabelPos) > 1
                warning('More than one match for specified reference channel; First one selected. This may cause eratic behavior. If the 2 channels are identical, delete one of them.');
            end
            if ~isempty( newLabelPos )
                icachansind(i) = newLabelPos(1);
                rminds(icachansind(i) == rminds) = [];
            else
                icachansind(i) = [];
            end
        end
        newICAchaninds(rminds,:) = [];
        EEG.icawinv = newICAchaninds;
        
        EEG.icachansind = icachansind;
        if length(EEG.icachansind) ~= size(EEG.icawinv,1)
            warning('Wrong channel indices, removing ICA decomposition');
            EEG.icaweights = [];
            EEG.icasphere  = [];
        else
            EEG.icaweights = pinv(EEG.icawinv);
            EEG.icasphere  = eye(length(icachansind));
        end;    
    end
    EEG = eeg_checkset(EEG);
end

% generate the output command
% ---------------------------
com = sprintf('EEG = pop_reref( EEG, %s);', vararg2str({ref, options{:}}));
