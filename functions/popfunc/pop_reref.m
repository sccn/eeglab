% pop_reref() - Convert an EEG dataset to average reference or to a
%               new common reference.
%
% Usage:
%       >> EEGOUT = pop_reref( EEG ); % pop up interactive window
%       >> EEGOUT = pop_reref( EEG, ref, 'key', 'val' ...);
%
% Graphical interface:
%   "Compute average reference" - [edit box] Checking this box is the same 
%                 as giving an empty value to the 'ref' command line parameter.
%   "Re-reference data to channel number" - [checkbox] Checking this option
%                 automatically uncheck the first checkbox and allow to enter a
%                 value in the edit box on the right of the checkbox. It does not
%                 correspond to any command line input.
%   "Re-reference data to channel number" - [edit box] Enter the index of the 
%                 electrode for re-referencing in this box. Same as using the 
%                 'ref' command line input.
%   "Include old reference channel" - [checkbox] When re-referencing the data,
%                 checking this checkbox allow to generate data for the old 
%                 reference channel. The location for this channel might not
%                 have been defined and can be specified using the following 
%                 textbox. Using this option is the same as setting the 'method'
%                 option to 'withref'.
%   "Include old reference channel" - [edit boxes] use these edit boxes to specify 
%                 the old reference channel location. Same as using the 'refloc' 
%                 optional input.  (note that if the channel location structure 
%                 contains one more channel that the data, the last channel is 
%                 considered to be the reference and the location for this channel
%                 does not have to be specified here). 
%
% Inputs:
%   EEG         - input dataset
%   ref         - reference: [] = average reference
%                            X  = new reference electrode number
%
% Optional inputs:
%   'method'    - ['standard'|'withref'] can be either 'standard' or 'withref' 
%                 to recompute the old reference potential. See also reref().
%   'refloc'    - [cell array] old common reference name polar location 
%                 (can also be included as the last channel of the EEG.chanlocs 
%                 struture). i.e. { 'cz' 0 0 }
%
% Inputs:
%   EEGOUT      - output dataset
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
    
    enable1 = fastif(strcmp(EEG.ref, 'averef'), 'on', 'off');
    enable2 = fastif(strcmp(EEG.ref, 'averef'), 'off', 'on');
    % build gui
	% ---------
    geometry = { [1] [1.8 1] [1.8 1] [1] [3 1 1 1] [3 1 1 1] };
    uilist = { { 'style' 'text' 'string' ['Data reference state is: ' EEG.ref] } ...
               { 'style' 'checkbox' 'tag' 'ave' 'value' fastif(strcmp(EEG.ref, 'averef'), 0, 1) 'string' 'Compute average reference' ...
                        'enable' fastif(strcmp(EEG.ref, 'averef'), 'off', 'on')  'callback' ...
                 [ 'set(findobj(''parent'', gcbf, ''tag'', ''reref''), ''value'', ~get(gcbo, ''value''));' ...
                   'set(findobj(''parent'', gcbf, ''tag'', ''rerefstr''), ''enable'', fastif(get(gcbo, ''value''), ''off'', ''on''));' ] } ...
               { } ...
               { 'style' 'checkbox' 'tag' 'reref' 'value' fastif(strcmp(EEG.ref, 'averef'), 1, 0) 'string' 'Re-reference data to channel number(s):' ...
                 'callback' [ 'set(findobj(''parent'', gcbf, ''tag'', ''ave''), ''value'', ~get(gcbo, ''value''));'  ...
                 'set(findobj(''parent'', gcbf, ''tag'', ''rerefstr''), ''enable'', fastif(get(gcbo, ''value''), ''on'', ''off''));' ] } ...
               { 'style' 'edit' 'tag' 'rerefstr' 'string' '' 'enable' enable1 } ...
               { } ...
               { } ...
               { 'style' 'text' 'tag' 'oldref' 'enable' 'off' 'string' 'Label' } ...
               { 'style' 'text' 'tag' 'oldref' 'enable' 'off' 'string' 'Theta' } ...
               { 'style' 'text' 'tag' 'oldref' 'enable' 'off' 'string' 'Radius' } ...
               { 'style' 'checkbox' 'enable' enable2 'string' 'Include old reference channel' 'callback' ...
                 'set(findobj(''parent'', gcbf, ''tag'', ''oldref''), ''enable'', fastif(get(gcbo, ''value''), ''on'', ''off''));' } ...
               { 'style' 'edit' 'tag' 'oldref' 'enable' 'off' 'string' '' } ...
               { 'style' 'edit' 'tag' 'oldref' 'enable' 'off' 'string' '' } ...
               { 'style' 'edit' 'tag' 'oldref' 'enable' 'off' 'string' '' } ...
             };
    if length(EEG.chanlocs) == EEG.nbchan+1
        geometry = { geometry{1:end-2}  [1] };
        uilist   = { uilist{1:end-8} ...
                     { 'style' 'checkbox' ...
                       'string' 'Include old reference channel (use location from the electrode location structure)' } };
    end;
    
    result = inputgui(geometry, uilist, 'pophelp(''pop_reref'')', 'pop_reref - average reference or re-reference data');

    % decode inputs
    % -------------
    if isempty(result), return; end;
    if ~isempty(result{3}), ref = eval([ '[' result{3} ']' ] );
    else                    ref = [];
    end;
    if result{1}, ref = []; end;
    options = { };
    if length(result) > 3 & result{4}, options = { options{:} 'method' 'withref' }; end;
    if length(result) > 4 & ~isempty(result{5}), 
        options = { options{:} 'refloc' { result{5} eval(result{6}) eval(result{7}) } };
    end;
else 
    options = varargin;
end;
optionscall = options;

% warning the user
% -----------------
withref = 0;
for index = 1:length(options)
    if isstr(options{index}) & strcmpi(options{index}, 'withref');
        withref = 1;
    end;
end;
if nargin < 2
    if length(ref) > 1 | withref == 0
        res = questdlg2(strvcat('Using multiple references, or using a single reference', ...
                                'without including the old reference reduces the dimensionality', ...
                                'of the data and prevents from re-referencing accuratelly later on.', ...
                                'Do you want to continue ?'), 'warning', 'Cancel', 'Yes', 'yes');
        if strcmpi(res, 'Cancel'), return; end;
    end;
end;

% include channel location file
% -----------------------------
if ~isempty(EEG.chanlocs)
    optionscall = { optionscall{:} 'elocs' EEG.chanlocs }; 
end;    
if isfield(EEG, 'ref')
    optionscall = { optionscall{:} 'refstate' EEG.ref }; 
end;

% include ICA or not
% ------------------
if ~isempty(EEG.icaweights)
    optionscall = { optionscall{:} 'icaweight' EEG.icaweights*EEG.icasphere }; 
    [EEG.data EEG.chanlocs EEG.icaweights EEG.icasphere] = reref(EEG.data, ref, optionscall{:});
else 
    [EEG.data EEG.chanlocs ] = reref(EEG.data, ref, optionscall{:});
end;

% add a tag in the dataset and clear some fields
% ----------------------------------------------
if isempty(ref)
    if strcmp(EEG.ref, 'averef')
        EEG.ref = 'common';
    else
        EEG.ref = 'averef';
    end;
else 
    EEG.ref = 'common';
end;
EEG.icaact  = [];
EEG.icawinv = [];
%EEG = eeg_checkset(EEG);

% generate the output command
% ---------------------------
if isempty( options )
    com = sprintf('%s = pop_reref( %s, [%s]);', inputname(1), inputname(1), num2str(ref));
else 
    com = sprintf('%s = pop_reref( %s, [%s], %s);', inputname(1), inputname(1), num2str(ref), vararg2str(options));
end;    
return;
