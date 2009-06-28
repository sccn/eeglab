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

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.33  2008/04/16 17:52:46  arno
% Additional entry to exclude reference from rereferencing
%
% Revision 1.32  2007/02/16 18:59:24  scott
% clarified help msg  -sm
%
% Revision 1.31  2006/05/04 10:07:07  arno
% same
%
% Revision 1.29  2005/01/24 19:30:37  arno
% remove field for Keun re-refencing
% ,.
%
% Revision 1.28  2004/05/14 23:58:23  arno
% operator precedence
%
% Revision 1.27  2004/01/30 23:01:52  arno
% update channel position only if refloc is set
%
% Revision 1.26  2003/11/05 16:24:17  arno
% homogenous -> homogeneous
%
% Revision 1.25  2003/10/14 17:03:08  arno
% vararg2str for optional arguments
%
% Revision 1.24  2003/09/08 22:54:12  arno
% typo for refstate
%
% Revision 1.23  2003/07/31 17:10:05  arno
% conversion to 3-D for the last channel
%
% Revision 1.22  2003/07/31 17:03:31  arno
% empty last channel
%
% Revision 1.21  2003/07/30 18:03:59  arno
% allowing empty channel location
%
% Revision 1.20  2003/07/28 17:53:39  arno
% channel ref index
%
% Revision 1.19  2003/07/28 16:46:49  arno
% remove redundancy
%
% Revision 1.18  2003/07/28 16:42:03  arno
% text with scott
%
% Revision 1.17  2003/07/27 01:19:21  arno
% debuging GUI call
%
% Revision 1.16  2003/07/27 01:09:30  arno
% debuging
%
% Revision 1.15  2003/07/27 00:44:32  arno
% brand new function for rereferencing
%
% Revision 1.14  2003/07/25 23:49:54  arno
% change interface, warning messages ...
%
% Revision 1.13  2003/07/25 00:11:56  arno
% allowing multiple references
%
% Revision 1.12  2003/07/25 00:10:51  arno
% correct typo
%
% Revision 1.11  2003/07/02 01:06:00  arno
% remove debug msg
%
% Revision 1.10  2003/07/02 01:05:34  arno
% debug msg
%
% Revision 1.9  2003/06/11 00:06:16  arno
% header message
%
% Revision 1.8  2003/02/17 02:53:22  arno
% reformating text for new functionality in help2html
%
% Revision 1.7  2003/02/16 23:30:20  arno
% adding GUI info
%
% Revision 1.6  2002/11/13 23:13:51  arno
% gui mutual exclusion problem
%
% Revision 1.5  2002/11/13 23:04:40  arno
% removing extra electrode in average ref
%
% Revision 1.4  2002/11/13 20:29:44  arno
% debugging
%
% Revision 1.3  2002/11/13 19:22:23  arno
% averef field -> ref field
%
% Revision 1.2  2002/11/12 23:23:37  arno
% mode -> method keyword
%
% Revision 1.1  2002/11/12 19:08:34  arno
% Initial revision
%
% Revision 1.1  2002/04/05 17:32:13  arno
% Initial revision
%

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
    cb_chansel1 = '[tmp tmpval] = pop_chansel({EEG(1).chanlocs.labels}, ''withindex'', ''on''); set(findobj(gcbf, ''tag'', ''reref''   ), ''string'',tmpval); clear tmp tmpval';
    cb_chansel2 = '[tmp tmpval] = pop_chansel({EEG(1).chanlocs.labels}, ''withindex'', ''on''); set(findobj(gcbf, ''tag'', ''exclude'' ), ''string'',tmpval); clear tmp tmpval';
    cb_chansel3 = [ 'if ~isfield(EEG(1).chaninfo, ''nodatchans''), ' ...
                    '   warndlg2(''There are no Reference channel defined, add it using the channel location editor'');' ...
                    'elseif isempty(EEG(1).chaninfo.nodatchans),' ...
                    '   warndlg2(''There are no Reference channel defined, add it using the channel location editor'');' ...
                    'else,' ...
                    '   [tmp tmpval] = pop_chansel({EEG(1).chaninfo.nodatchans.labels}, ''withindex'', ''on''); set(findobj(gcbf, ''tag'', ''refloc''  ), ''string'',tmpval); clear tmp tmpval;' ...
                    'end;' ];
    
    % find current reference (= reference most used)
    % ----------------------------------------------
    if isfield(EEG(1).chanlocs, 'ref')
        [curref tmp allinds] = unique( { EEG(1).chanlocs.ref });
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
             tmpallchans = lower({ EEG.chaninfo.nodatchans.labels });
             chanind = strmatch( lower(restag.refloc), tmpallchans, 'exact');
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
[EEG.data EEG.chanlocs inds ] = reref(EEG.data, ref, optionscall{:});
EEG.nbchan = size(EEG.data,1);
EEG = eeg_checkset(EEG);

% include ICA or not
% ------------------
if ~isempty(EEG.icaweights)
    g = struct(optionscall{:});
    if ~isfield(g, 'exclude'), g.exclude = []; end;
    if ~isfield(g, 'keepref'), g.keepref = 'off'; end;
    if ~isfield(g, 'refloc') , g.refloc  = []; end;
    
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
        EEG.icawinv = reref(EEG.icawinv, ref, optionscall{:});
        
        % get output channel indices
        % --------------------------
        chansout = 1:nchans;
        if ~isempty(ref) & strcmpi(g.keepref,'off')
            ref = sort(ref);
            for ind = length(ref):-1:1
                chansout(ref(ind)+1:end) = chansout(ref(ind)+1:end)-1;
                chansout(ref(ind)) = [];
            end;
        end;
        
        % convert channel indices in icachanlocs
        % --------------------------------------
        icachansind = EEG.icachansind;
        for i=length(icachansind):-1:1
            indchan = find( icachansind(i) == chansout );
            if ~isempty( indchan )
                icachansind(i) = indchan;
            else
                icachansind(i) = [];
            end;
        end;
        
        % add new channel if necessary
        if ~isempty(g.refloc)
            icachansind = [ icachansind size(EEG.data,1) ];
        end;
        
        EEG.icachansind = icachansind;
        if length(EEG.icachansind) ~= size(EEG.icawinv,1)
            warning('Wrong channel indices, removing ICA decomposition');
            dsafdsf
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
