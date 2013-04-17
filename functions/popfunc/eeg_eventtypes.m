% eeg_eventtypes()  - return a list of event or urevent types in a dataset and 
%                     the respective number of events of each type. Ouput event 
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
%           >> eeg_eventtypes(EEG);       % print numner of each event types
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
%           >> eeg_eventtypes(EEG);       % print numner of each event types
%
% Curently disabled:
%           >> eeg_eventtypes(EEG,'urevent');  % print hist. of urevent types
%           >> eeg_eventtypes(EEG,{'rt'});% print number of 'rt' events 
%           >> eeg_eventtypes(EEG,'urevent',{'rt','break'}); 
%                                         % print numbers of 'rt' and 'break' 
%                                         % type urevents 

% Copyright (C) 2004 Scott Makeig, SCCN/INC/UCSD, smakeig@ucsd.edu
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
end;

UREVENTS = 0; % flag returning infor for urevents instead of events
typelist = [];
if nargin>1
   if ischar(arg2)
       if strcmp(arg2,'urevent') | strcmp(arg2,'urevents')
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
