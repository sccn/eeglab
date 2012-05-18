% eeg_epochformat() - Convert the epoch information of a dataset from struct
%                     to array or vice versa.
%
% Usage: >> [epochsout fields] = eeg_epochformat( epochs, 'format', fields, events );
%
% Input:
%   epochs - epoch numerical or cell array or epoch structure
%   format - ['struct'|'array'] convert epoch array to structure and epoch
%            structure to array.
%   fields - [optional] cell array of strings containing the names of
%            the epoch struct fields. If this field is empty, it uses the
%            following list for the names of the fields { 'var1' 'var2' ... }.
%            For structure conversion, this field helps export a given
%            event type. If this field is left empty, the time locking
%            event for each epoch is exported.
%   events - numerical array of event indices associated with each epoch.
%            For array conversion, this field is ignored.
%
% Outputs:
%   epochsout - output epoch array or structure
%   fields    - output cell array with the name of the fields
%
% Epoch format:
%   struct - Epoch information is organised as an array of structs
%   array  - Epoch information is organised as an 2-d array of numbers, 
%            each column representing a user-defined variable (the 
%            order of the variable is a function of its order in the 
%            struct format).
%
% Note: 1) The epoch structure is defined only for epoched data.
%       2) The epoch 'struct' format is more comprehensible.
%          For instance, to see all the properties of epoch i,
%          type >> EEG.epoch(i)
%          Unfortunately, structures are awkward for expert users to deal
%          with from the command line (Ex: To get an array of 'var1' values,
%           >> celltomat({EEG.epoch(:).var1})')
%          In array format, asuming 'var1' is the first variable
%          declared, the same information is obtained by
%           >> EEG.epoch(:,1)
%       3) This function automatically updates the 'epochfields'
%          cell array depending on the format.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 12 Feb 2002
%
% See also: eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) from eeg_eventformat.m, 
% Arnaud Delorme, CNL / Salk Institute, 12 Feb 2002, arno@salk.edu
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

% $Log: eeg_epochformat.m,v $
% Revision 1.4  2005/05/24 16:57:09  arno
% cell2mat
%
% Revision 1.3  2003/07/20 19:32:20  scott
% typos
%
% Revision 1.2  2002/04/21 01:10:35  scott
% *** empty log message ***
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 03/13/02 added field arrays options -ad

function [ epoch, fields, epocheventout] = eeg_epochformat( epoch, format, fields, epochevent);

if nargin < 2
   help eeg_epochformat;
   return;
end;	

if nargin < 3
    fields = {};
end;
epocheventout = [];
    
switch format
case 'struct'
   if ~isempty(epoch) & ~isstruct(epoch)

      fields = getnewfields( fields, size(epoch,2) - length(fields));   
       
      % generate the structure
      % ----------------------
      command = 'epoch = struct(';
      for index = 1:length(fields)
          if iscell(epoch)
              command = [ command '''' fields{index} ''', epoch(:,' num2str(index) ')'',' ];
          else
              command = [ command '''' fields{index} ''', mattocell( epoch(:,' num2str(index) ')'',' ...
                            '[1], ones(1,size(epoch,1))),' ];
          end;                  
      end;
      eval( [command(1:end-1) ');' ] );

      if exist('epochevent') == 1
         for index = 1:size(epoch,2)
             if iscell(epochevent)
                 epoch(index).event = epochevent{index};
             else
                 epoch(index).event = epochevent(index);
             end;
         end;
      end;
   end

case 'array'
    if isstruct(epoch)
        
        % note that the STUDY std_maketrialinfo also gets the epoch info for the
        % time locking event
        
        selectedType = fields;
        if iscell(fields) && ~isempty(fields), selectedType = fields{1}; end;
  	    fields = fieldnames( epoch );
        
        eval( [ 'values = { epoch.' fields{1} ' };' ]);
        
        if any(cellfun(@length, values) > 1)
            if ~isfield(epoch, 'eventlatency')
                error('eventlatency field not present in data epochs');
            end;
            
            if isempty(selectedType)
                % find indices of time locking events
                for index = 1:length(epoch)
                    epochlat = [ epoch(index).eventlatency{:} ];
                    tmpevent = find( abs(epochlat) < 0.02 );
                    if isempty(tmpevent)
                        error('time locking event missing, cannot convert to array');
                    end;
                    epochSubIndex(index) = tmpevent;
                end;
            else
                % find indices of specific event type (if several take the
                % first one
                for index = 1:length(epoch)
                    epochtype = epoch(index).eventtype;
                    tmpeventind = strmatch( selectedType, epochtype );
                    if length(tmpeventind) > 1
                        fprintf('Warning: epoch %d has several events of "type" %s, taking the fist one\n', index, selectedType);
                    end;
                    if isempty(tmpeventind)
                         epochSubIndex(index) = NaN;
                    else epochSubIndex(index) = tmpeventind(1);
                    end;
                end;
            end;
        else
            epochSubIndex = ones(1, length(epoch));
        end;
        
        % copy values to array
        tmp = cell( length(epoch), length( fields ));
        for index = 1:length( fields )
            for trial = 1:length(epoch)
                tmpval = getfield(epoch, {trial}, fields{index});
                if isnan(epochSubIndex(trial))
                     tmp(trial, index) = { NaN };
                elseif iscell(tmpval)
                     tmp(trial, index) = tmpval(epochSubIndex(trial));
                elseif ~ischar(tmpval)
                     tmp(trial, index) = { tmpval(epochSubIndex(trial)) };
                else tmp(trial, index) = { tmpval };
                end;
            end;
        end;
        epoch = tmp;
    end;
    otherwise, error('unrecognised format');   
end;

return;

% create new field names
% ----------------------
function epochfield = getnewfields( epochfield, nbfields )
   count = 1;
   if nbfields > 0
        while nbfields > 0
            if isempty( strmatch([ 'var' int2str(count) ], epochfield ) )
               epochfield =  { epochfield{:} [ 'var' int2str(count) ] };
               nbfields = nbfields-1;
            else    count = count+1;
            end;
        end;                        
   else
        epochfield = epochfield(1:end+nbfields);
   end;     
return;

