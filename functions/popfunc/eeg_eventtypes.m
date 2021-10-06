% eeg_eventtypes()  - return a list of event or urevent types in a dataset and 
%                     the respective number of events of each type. Output event 
%                     types are sorted in reverse order of their number. If no 
%                     outputs, print this list on the commandline instead.
%
% Usage:
%        >> [types,numbers] = eeg_eventtypes(EEG);
% Inputs:
%        EEG        - EEGLAB dataset structure
% Outputs:
%        types      - cell array of event type strings
%        numbers    - vector giving the numbers of each event type in the data
%
% Example:
%           >> eeg_eventtypes(EEG);       % print number of each event types
%
% Author: Scott Makeig, SCCN/INC/UCSD, April 28, 2004-

%        >> [types,numbers] = eeg_eventtypes(EEG,types);
%        >> [types,numbers] = eeg_eventtypes(EEG,'urevents',types);
% Inputs:
%        EEG        - EEGLAB dataset structure
%        'urevents' - return event information for the EEG.urevent structure
%        types      - {cell array} of event types to return or print.
% Outputs:
%        types      - cell array of event type strings
%        numbers    - vector giving the numbers of each event type in the data
%
% Note:  Numeric (ur)event types are converted to strings, so, for example, 
%        types {13} and {'13'} are not distinguished.
%
% Example:
%           >> eeg_eventtypes(EEG);       % print number of each event types
%
% Currently disabled:
%           >> eeg_eventtypes(EEG,'urevent');  % print hist. of urevent types
%           >> eeg_eventtypes(EEG,{'rt'});% print number of 'rt' events 
%           >> eeg_eventtypes(EEG,'urevent',{'rt','break'}); 
%                                         % print numbers of 'rt' and 'break' 
%                                         % type urevents 

% Copyright (C) 2004 Scott Makeig, SCCN/INC/UCSD, smakeig@ucsd.edu
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
%
%
% event types can be numbers, Stefan Debener, 05/12/2006, 
% added 2nd and 3rd args, sorted outputs by number, Scott Makeig, 09/12/06

function [types,numbers] = eeg_eventtypes(EEG,arg2,arg3)

if nargin< 1
   help eeg_eventtypes
   return
end
if ~isstruct(EEG)
   error('EEG argument must be a dataset structure')
end

if ~isfield(EEG,'event')
   error('EEG.event field not found');
end
if nargin > 1
    error('Multiple input arguments are currently disabled');
end

UREVENTS = 0; % flag returning info for urevents instead of events
typelist = [];
if nargin>1
   if ischar(arg2)
       if strcmp(arg2,'urevent') || strcmp(arg2,'urevents')
            UREVENTS = 1; % change flag
       else
            error('second argument string not understood')
       end
       if nargin>2
            if iscell(arg3)
                  typelist = arg3;
            end
       end
   elseif iscell(arg2)
       typelist = arg2;
   end
end
if ~isempty(typelist)     % cast to cell array of strings
   for k=1:length(typelist)
       if isnumeric(typelist{k})
          typelist{k} = num2str(typelist{k});
       end
   end
end
       
if ~UREVENTS
   nevents = length(EEG.event);
   alltypes = cell(nevents,1);
   for k=1:nevents
       if isnumeric(EEG.event(k).type)
          alltypes{k} = num2str(EEG.event(k).type);
       else
          alltypes{k} = EEG.event(k).type;
       end
   end
else
   nevents = length(EEG.urevent);
   alltypes = cell(nevents,1);
   for k=1:nevents
       if isnumeric(EEG.urevent(k).type)
          alltypes{k} = num2str(EEG.urevent(k).type);
       else
          alltypes{k} = EEG.urevent(k).type;
       end
   end
end

[types i j] = unique_bc(alltypes);

istypes = 1:length(types);
notistypes = [];
if ~isempty(typelist)
  notistypes = ismember_bc(typelist,types);
  istypes = ismember_bc(types,typelist(find(notistypes==1))); % types in typelist?
  notistypes = typelist(find(notistypes==0));
end

types(~istypes) = []; % restrict types to typelist

ntypes = length(types);
numbers = zeros(ntypes + length(notistypes),1);
for k=1:ntypes
  numbers(k) = length(find(j==k));
end
types = [types(:); notistypes(:)]; % concatenate the types not found
ntypes = length(types);
for j = 1:length(notistypes)
   numbers(k+j) = 0;
end

% sort types in reverse order of event numbers
[numbers nsort] = sort(numbers);
numbers = numbers(end:-1:1);
types = types(nsort(end:-1:1)); 

%
% print output on commandline %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if nargout < 1 
  fprintf('\n');
  if UREVENTS
    fprintf('EEG urevent types:\n\n')
  else
    fprintf('EEG event types:\n\n')
  end

  maxx = 0;
  for k=1:ntypes
    x = length(types{k});
    if x > maxx
        maxx = x; % find max type name length
    end
  end
  for k=1:ntypes
    fprintf('  %s',types{k});
    for j=length(types{k})+1:maxx+4
       fprintf(' ');
    end
    fprintf('%d\n',numbers(k));
 end
 fprintf('\n');
 clear types % no return variables
end
