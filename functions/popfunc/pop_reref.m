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
%   "Retain old reference channels in data" - [checkbox] When re-referencing the 
%                 data, checking this checkbox includes the data for the 
%                 previous reference channel.
%   "Exclude channel indices (EMG, EOG)" - [edit box] exclude the given
%                 channel indices from rereferencing.
%   "Add current reference channel back to the data" - [edit box] When 
%                 re-referencing the data, checking this checkbox
%                 reconstitutes the data for the previous reference
%                 channel. If the location for this channel  was not 
%                 defined, it can be specified using the text box below.
% Inputs:
%   EEG         - input dataset
%   ref         - reference: []            = convert to average reference
%                            [int vector]  = new reference electrode number(s)
% Optional inputs:
%   'exclude'   - [integer array] List of channels to exclude.
%   'keepref'   - ['on'|'off'] keep the reference channel.
%   'refloc'    - [structure] Previous reference channel structure.
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

function [EEG, com] = pop_reref( EEG, ref, varargin);

com = '';
if nargin < 1
   help pop_reref;
   return;
end;   
if isempty(EEG.data)
    error('Pop_reref: cannot process empty data');
end;

% gui inputs
% ----------
orichanlocs = EEG.chanlocs;
orinbchan   = EEG.nbchan;
if nargin < 2
    
    % find initial reference
    % ----------------------
    if length(EEG.chanlocs) == EEG.nbchan+1
        includeref = 1;
    end;
    
    geometry = { [1] [1] [1.8 1 0.3] [1] [1] [1.8 1 0.3] [1.8 1 0.3] };
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
    cb_chansel1 = 'tmpchanlocs = EEG(1).chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels}, ''withindex'', ''on''); set(findobj(gcbf, ''tag'', ''reref''   ), ''string'',tmpval); clear tmpchanlocs tmp tmpval';
    cb_chansel2 = 'tmpchanlocs = EEG(1).chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels}, ''withindex'', ''on''); set(findobj(gcbf, ''tag'', ''exclude'' ), ''string'',tmpval); clear tmpchanlocs tmp tmpval';
    cb_chansel3 = [ 'if ~isfield(EEG(1).chaninfo, ''nodatchans''), ' ...
                    '   warndlg2(''There are no Reference channel defined, add it using the channel location editor'');' ...
                    'elseif isempty(EEG(1).chaninfo.nodatchans),' ...
                    '   warndlg2(''There are no Reference channel defined, add it using the channel location editor'');' ...
                    'else,' ...
                    '   tmpchaninfo = EEG(1).chaninfo; [tmp tmpval] = pop_chansel({tmpchaninfo.nodatchans.labels}, ''withindex'', ''on''); set(findobj(gcbf, ''tag'', ''refloc''  ), ''string'',tmpval); clear tmpchanlocs tmp tmpval;' ...
                    'end;' ];
    if isempty(EEG.chanlocs), cb_chansel1 = ''; cb_chansel2 = ''; cb_chansel3 = ''; end;
    
    % find current reference (= reference most used)
    % ----------------------------------------------
    if isfield(EEG(1).chanlocs, 'ref')
        tmpchanlocs = EEG(1).chanlocs;
        [curref tmp allinds] = unique( { tmpchanlocs.ref });
        maxind = 1;
        for ind = unique(allinds)
            if length(find(allinds == ind)) > length(find(allinds == maxind))
                maxind = ind;
            end;
        end;
        curref = curref{maxind};
        if isempty(curref), curref = 'unknown'; end;
    else curref = 'unknown';
    end;
    
    uilist = { { 'style' 'text' 'string' [ 'Current data reference state is: ' curref] } ...
               ...
               { 'style' 'checkbox' 'tag' 'ave'   'value' 1 'string' 'Compute average reference' 'callback' cb_averef } ...
               ...
               { 'style' 'checkbox' 'tag' 'rerefstr' 'value' 0 'string' 'Re-reference data to channel(s):' 'callback'  cb_ref } ...
               { 'style' 'edit' 'tag' 'reref' 'string' '' 'enable' 'off' } ...
               { 'style' 'pushbutton' 'string' '...' 'callback' cb_chansel1 'enable' 'off' 'tag' 'refbr' } ...
               ...
               {} ...
               ...
               { 'style' 'checkbox' 'value' 0 'enable' 'off' 'tag' 'keepref' 'string' 'Retain old reference channels in data' } ...
               ...
               { 'style' 'text' 'string' 'Exclude channel indices (EMG, EOG)' } ...
               { 'style' 'edit' 'tag' 'exclude' 'string' '' } ...
               { 'style' 'pushbutton' 'string' '...' 'callback' cb_chansel2 } ...
               ...
               { 'style' 'text' 'tag' 'reflocstr' 'string' 'Add current reference channel back to the data' } ...
               { 'style' 'edit' 'tag' 'refloc' 'string' '' } ...
               { 'style' 'pushbutton' 'string' '...' 'callback' cb_chansel3 } };
    
    [result tmp tmp2 restag] = inputgui(geometry, uilist, 'pophelp(''pop_reref'')', 'pop_reref - average reference or re-reference data');
    if isempty(result), return; end;

    % decode inputs
    % -------------
    options = {};
    if ~isempty(restag.refloc),
         try
             tmpchaninfo = EEG.chaninfo;
             tmpallchans = lower({ tmpchaninfo.nodatchans.labels });
             allelecs = parsetxt(lower(restag.refloc));
             chanind  = [];
             for iElec = 1:length(allelecs)
                 chanind = [chanind strmatch( allelecs{iElec}, tmpallchans, 'exact') ];
             end;
             options = { options{:} 'refloc' EEG.chaninfo.nodatchans(chanind) }; 
         catch, disp('Error with old reference: ignoring it');
         end;
    end;
    if ~isempty(restag.exclude), options = { options{:} 'exclude' eeg_chaninds(EEG, restag.exclude) }; end;
    if restag.keepref,           options = { options{:} 'keepref' 'on' }; end;
    if restag.ave,               ref = []; end;
    if restag.rerefstr           
        if isempty(restag.reref)
            warndlg2('Abording: you must enter one or more reference channels'); 
            return;
        else
            ref = eeg_chaninds(EEG, restag.reref); 
        end;
    end;
else
    options = varargin;
end;
optionscall = options;

% include channel location file
% -----------------------------
if ~isempty(EEG.chanlocs)
    optionscall = { optionscall{:} 'elocs' EEG.chanlocs }; 
end;    

nchans = EEG.nbchan;
fprintf('Re-referencing data\n');
oldchanlocs = EEG.chanlocs;
[EEG.data EEG.chanlocs refchan ] = reref(EEG.data, ref, optionscall{:});
g = struct(optionscall{:});
if ~isfield(g, 'exclude'), g.exclude = []; end;
if ~isfield(g, 'keepref'), g.keepref = 'off'; end;
if ~isfield(g, 'refloc') , g.refloc  = []; end;

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
        for ind = 1:length(allf)
            EEG.chaninfo.nodatchans = setfield(EEG.chaninfo.nodatchans, { n }, ...
                allf{ind}, getfield(refchan, allf{ind}));
        end;
    end;
end;
if ~isempty(g.refloc)
    allinds = [];
    tmpchaninfo = EEG.chaninfo;
    for iElec = 1:length(g.refloc)
         allinds = [allinds strmatch( g.refloc(iElec).labels, { tmpchaninfo.nodatchans.labels }) ];
    end;
    EEG.chaninfo.nodatchans(allinds) = [];
end;
    
% legacy EEG.ref field
% --------------------
if isfield(EEG, 'ref')
    if strcmpi(EEG.ref, 'common') && isempty(ref)
        EEG.ref = 'averef';
    elseif strcmpi(EEG.ref, 'averef') && ~isempty(ref)
        EEG.ref = 'common';
    end;
end;

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
        end;
        
        newICAchaninds = zeros(orinbchan, size(EEG.icawinv,2));
        newICAchaninds(EEG.icachansind,:) = EEG.icawinv;
        
        [newICAchaninds newchanlocs] = reref(newICAchaninds, ref, optionscall{:});
        
        % convert channel indices in icachanlocs (uses channel labels)
        % ------------------------------------------------------------
        icachansind = EEG.icachansind;
        rminds      = [1:size(newICAchaninds,1)];
        for i=length(icachansind):-1:1
            oldLabel    = orichanlocs(icachansind(i)).labels;
            newLabelPos = strmatch(oldLabel, { newchanlocs.labels }, 'exact');
            
            if ~isempty( newLabelPos )
                icachansind(i) = newLabelPos;
                rminds(find(icachansind(i) == rminds)) = [];
            else
                icachansind(i) = [];
            end;
        end;
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
    end;
    EEG = eeg_checkset(EEG);
end;

% generate the output command
% ---------------------------
com = sprintf('%s = pop_reref( %s, %s);', inputname(1), inputname(1), vararg2str({ref, options{:}}));
