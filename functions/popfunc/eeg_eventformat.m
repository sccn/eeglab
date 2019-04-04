% eeg_eventformat() - Convert the event information of a dataset from struct
%                 to array or vice versa.
%
% Usage: >> [eventout fields] = eeg_eventformat( event, 'format', fields );
%
% Inputs:
%   event  - event array or structure
%   format - ['struct'|'array'] see below
%   fields - [optional] cell array of strings containing the names of
%            the event struct fields. If this field is empty, it uses 
%            the following list for
%            the names of the fields { 'type' 'latency' 'var1' ...
%            'var2' ... }.
% Output:
%   eventout  - output event array or structure
%   fields    - output cell array with the name of the fields
%
% Event formats:
%   struct - Events are organised as an array of structs with at
%            least two fields ('type' and 'latency')
%            (Ex: reaction_time may be type 1).
%   array  - events are organized as an array, the first column
%            representing the type, the second the latency and the
%            other ones user-defined variables.
%
% Note: 1) The event structure is defined only for continuous data
%          or epoched data derived from continuous data.
%       2) The event 'struct' format is more comprehensible.
%          For instance, to see all the properties of event 7,
%          type >> EEG.event(7)
%          Unfortunately, structures are awkward for expert users to deal
%          with from the command line (Ex: To get an array of latencies,
%           >> [ EEG.event(:).latency ] )
%          In array format, the same information is obtained by typing
%           >> EEG.event(:,2)
%       3) This function automatically updates the 'eventfield'
%          cell array depending on the format.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 27 Jan 2002
%
% See also: eeglab(), pop_selectevent(), pop_importevent()

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 27 Jan 2002, arno@salk.edu
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

% 2/06/02 modifed header - sm & ad
% 2/08/02 add field input - ad
% 2/12/02 reprogrammed function using epochformat.m - ad

function [event, eventfield] = eeg_eventformat(event, format, fields);

if nargin < 2
   help eeg_eventformat;
   return;
end;	

if exist('fields') ~= 1, fields = { 'type', 'latency' }; end

[event eventfield] = eeg_epochformat( event, format, fields);

