% pop_importpres() - append Presentation event file information into an EEGLAB dataset
%                    The Presentation stimulus presentation program outputs an ascii
%                    log file. This function merges existing EEG dataset events with
%                    additional field information (fields) about those events contained 
%                    in the logfile. 
% Usage:
%   >> EEGOUT = pop_importpres( EEGIN, filename );
%   >> EEGOUT = pop_importpres( EEGIN, filename, typefield, ...
%                                  latfield, durfield, align, 'key', 'val', ... );
% Inputs:
%   EEGIN          - input dataset
%   logfilename    - Presentation logfile name
%
%   typefield      - [string] type fieldname {default: 'code'}
%   latfield       - [string] latency fieldname {default: 'time'}
%   durfield       - [string] duration fieldname {default: 'none'}
%   align          - [integer] alignment with pre-existing events
%                    See    >> help pop_importevent
%   'key','val'    - This function calls pop_importevent(). These are
%                    optional arguments for this function (for event 
%                    alignment for instance).
% Outputs:
%   EEGOUT         - data structure with added Presentation logfile information
%
% Note: If there are pre-existing events in the input dataset,
%       this function will recalculate the latencies of the events
%       in the Presentation file, so that they match those
%       of the pre-existing events.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 15 March 2002
%
% Note: This function is backward compatible with its early versions
%       (before the input argument 'durfield' was introduced). 
%       It can read the 'align' value as its 5th (not 6th) paramater. 
%
% See also: eeglab(), pop_importevent()

% Copyright (C) 13 March 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [EEG, command] = pop_importpres(EEG, filename, typefield, latfield, durfield, align, varargin); 

command = '';

if nargin < 1 
    help pop_importpres;
    return
end

% decode input (and backward compatibility)
% -----------------------------------------
if nargin < 5
    durfield = '';
end
if nargin >= 5 && ~ischar(durfield)
    if nargin >= 6
        varargin = { align varargin{:} };
    end
    align = durfield;
    durfield = '';
else
    if nargin < 6
        align = 0;
    end
end

if nargin < 2 
	% ask user
	[filename, filepath] = uigetfile('*.log;*.LOG', 'Choose a Presentation file -- pop_importpres()'); 
    drawnow;
	if filename == 0 return; end
	filename = [filepath filename];
end

fields = loadtxt(filename, 'delim', 9, 'skipline', -2, 'nlines', 1, 'verbose', 'off');

% finding fields
% --------------
if nargin > 1
    if nargin < 3
        typefield = 'code'; % backward compatibility
        latfield  = 'time';
    end
    indtype  = strmatch(lower(typefield), lower(fields));
    indlat   = strmatch(lower(latfield) , lower(fields));
    if ~isempty(durfield)
         inddur   = strmatch(lower(durfield) , lower(fields));
    else inddur = 0;
    end
else
    indtype1   = strmatch('event type', lower(fields));
    indtype2   = strmatch('code', lower(fields));
    indtype    = [ indtype1 indtype2 1];
    indlatency = strmatch('time', lower(fields), 'exact');
    indlatency = [ indlatency 1 ];
    uilist = { { 'style' 'text' 'string' [ 'File field containing event types' 10 '' ] } ...
               { 'style' 'list' 'string' strvcat(fields) 'value' indtype(1)  'listboxtop' indtype(1)} ...
               { 'style' 'text' 'string' [ 'File field containing event latencies' 10 '' ] } ...
               { 'style' 'list' 'string' strvcat(fields) 'value' indlatency(1) 'listboxtop' indlatency(1) } ...
               { 'style' 'text' 'string' [ 'File field containing event durations' 10 '' ] } ...
               { 'style' 'list' 'string' strvcat({ 'None' fields{:} }) 'value' 1 'listboxtop' 1 } ...
               { } { 'style' 'text' 'string' 'Note: scroll lists then click to select field' } };
    uigeom = { [2 1] [2 1] [2 1] 1 1 };
    result = inputgui(uigeom, uilist, 'pophelp(''pop_importpres'')', 'Import presentation file - pop_importpres()', ...
                      [], 'normal', [2 2 2 1 1]);
    if isempty(result), return; end
    
    indtype = result{1};
    indlat  = result{2};
    inddur  = result{3}-1;
    typefield = fields{indtype};
    latfield  = fields{indlat};
    if inddur ~= 0
         durfield  = fields{inddur};
    else durfield  = '';
    end
end
if isempty(indtype)
    error(['Could not detect field ''' typefield ''', try importing the file as ASCII (use delimiter=9 (tab))']);
end
if isempty(indlat)
    error(['Could not detect field ''' latfield ''', try importing the file as ASCII (use delimiter=9 (tab))']);
end
disp(['Replacing field ''' typefield ''' by ''type'' for EEGLAB compatibility']);
disp(['Replacing field ''' latfield  ''' by ''latency'' for EEGLAB compatibility']);
fields{indtype} = 'type';
fields{indlat}  = 'latency';
if inddur ~= 0
    fields{inddur}  = 'duration';
end

% check inputs


% regularizing field names
% ------------------------
for index = 1:length(fields)
    indspace = find(fields{index} == ' ');
    fields{index}(indspace) = '_';
    indparen = find(fields{index} == ')');
    if ~isempty(indparen) && indparen == length(fields{index})
        % remove text for parenthesis
        indparen = find(fields{index} == '(');
        if indparen ~= 1
            disp([ 'Renaming ''' fields{index} ''' to ''' fields{index}(1:indparen-1) ''' for Matlab compatibility' ]);
            fields{index} = fields{index}(1:indparen-1);
        else
            fields{index}(indspace) = '_';
        end
    else
        fields{index}(indspace) = '_';
        indparen = find(fields{index} == '(');
        fields{index}(indspace) = '_';
    end
end

% find if uncertainty is duplicated
% ---------------------------------
induncert  = strmatch('uncertainty', lower(fields), 'exact');
if length(induncert) > 1
    fields{induncert(2)}= 'Uncertainty2';
    disp('Renaming second ''Uncertainty'' field');
end

% import file
% -----------
if isempty(EEG.event), align = NaN; end

%EEG = pop_importevent(EEG, 'append', 'no', 'event', filename, 'timeunit', 1E-4, 'skipline', -3, ...
%                           'delim', 9, 'align', align, 'fields', fields, varargin{:});
EEG = pop_importevent(EEG, 'event', filename, 'timeunit', 1E-4, 'skipline', -3, ...
                           'delim', 9, 'align', align, 'fields', fields, varargin{:});

command = sprintf('EEG = pop_importpres(EEG, %s);', vararg2str({ filename typefield latfield durfield align })); 

return;
