% pop_VisEd() - Collect variables for visual editing.
%
% Usage: 
%   >>  EEG = pop_VisEd( EEG, ChanIndex, EventType);
%
% Inputs:
%   ChanIndex   - EEG channels to display in eegplot figure window while editing events and identifying bad channels.
%   EventType   - Event types to display in eegplot figure window while editing events and identifying bad channels.
% 
% Optional Inputs:
%   quick_evtmk - When quick_evtmk is a string the crtl click command
%       bypasses the pop_fig_EditEvent GUI and simply adds a new event. The
%       string entered for the quick_evtmk input is used as the new event type
%       string. Default quick_evtmk = ''.
%   quick_evtrm - When quick_evtrm = 'on' the ctrl click command (when
%       selecting an existing event) removes the event without initiating the
%       pop_fig_EditEvent GUI. Default quick_evtrm = 'off'.
%
% Outputs:
%   EEG  - output dataset
%
% UI for selecting EEG channels and event types to be displayed in eegplot
% figure window while editing events and identifying bad channels.
%
% Calls function EEG=VisEd(EEG,ChanIndex,EventType);
%
% See also:
%   EEGLAB 

% Copyright (C) <2008>  <James Desjardins> Brock University
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

% Edit log:
%
% 2012 02 19; Optional inputs quick_evtmk and quick_evtrm added.

function [EEG,com]=pop_VisEd(EEG, ChanIndex, EventType)


% the command output is a hidden output that does not have to
% be described in the header
com = ''; % this initialization ensure that the function will return something
          % if the user press the cancel button            
          % display help if not enough arguments
% ------------------------------------

if nargin < 1
	help pop_VisEd;
	return;
end;	

DataTypeCell={'EEG'};
if ~isempty(EEG.icaweights);
    DataTypeCell={'EEG','ICA'};
    for i=1:length(EEG.icaweights(:,1));EEG.ic(i).labels=sprintf('%s%s','comp',num2str(i));end
end

if isempty(EEG.chanlocs);
    disp('Labelling channels by number.');
    for i=1:EEG.nbchan;
        EEG.chanlocs(i).labels=num2str(i);
    end
end

Num2StrEvCount=0;
for i=1:length(EEG.event);
    if isnumeric(EEG.event(i).type);
        Num2StrEvCount=Num2StrEvCount+1;
        EEG.event(i).type=num2str(EEG.event(i).type);
    end
end
if Num2StrEvCount>0;
    disp(sprintf('%s%s', num2str(Num2StrEvCount), 'event types converted from num2str'));
end


% pop up window
% -------------
if nargin < 3

    if ~isempty(EEG.event)
         tmpevent  = EEG.event;
         eventlist = vararg2str(unique({tmpevent.type}));
    else eventlist = '';
    end;
    results=inputgui( ...
    {[1] [1] [4 4 1] [4 4 1] [4 4 1] [1] [1] [1] [1]}, ...
    {...
        ... %1
        {'Style', 'text', 'string', 'Enter visual editing parameters.', 'FontWeight', 'bold'}, ...
        ... %2
        {}, ...
        ... %3
        {'Style', 'text', 'string', 'Data type to display:'}, ...
        {'Style', 'popup', 'string', DataTypeCell, 'tag', 'DataTypePop'... 
                  'callback', ['switch get(findobj(gcbf, ''tag'', ''DataTypePop''), ''value'');' ...
                               '    case 1;' ...
                               '        tmpchanlocs = EEG.chanlocs;' ...
                               '        set(findobj(gcbf, ''tag'', ''ChanLabelButton''), ''callback'',' ...
                               '            [''[ChanLabelIndex,ChanLabelStr,ChanLabelCell]=pop_chansel({tmpchanlocs.labels});' ...
                               '             set(findobj(gcbf, ''''tag'''', ''''ChanIndexEdit''''), ''''string'''', vararg2str(ChanLabelIndex))'']);' ...
                               '        set(findobj(gcbf, ''tag'', ''ChanIndexEdit''), ''string'', vararg2str(1:EEG.nbchan));' ...
                               '    case 2;' ...
                               '        set(findobj(gcbf, ''tag'', ''ChanLabelButton''), ''callback'',' ...
                               '            [''for i=1:length(EEG.icaweights(:,1));IC(i).labels=sprintf(''''%s%s'''',''''comp'''',num2str(i));end;' ...
                               '            [ChanLabelIndex,ChanLabelStr,ChanLabelCell]=pop_chansel({IC.labels});' ...
                               '             set(findobj(gcbf, ''''tag'''', ''''ChanIndexEdit''''), ''''string'''', vararg2str(ChanLabelIndex))'']);' ...
                               '        set(findobj(gcbf, ''tag'', ''ChanIndexEdit''), ''string'', vararg2str(1:length(EEG.icaweights(:,1))));' ...
                               'end; clear tmpchanlocs;']}, ...
        {}, ...
        ... %4
        {'Style', 'text', 'string', 'Channels to display in eegplot figure window:'}, ...
        {'Style', 'edit', 'string', vararg2str(1:EEG.nbchan),'tag', 'ChanIndexEdit'}, ...
        {'Style', 'pushbutton', 'string', '...', 'tag', 'ChanLabelButton',... 
                  'callback', ['tmpchanlocs = EEG.chanlocs; [ChanLabelIndex,ChanLabelStr,ChanLabelCell]=pop_chansel({tmpchanlocs.labels}); clear tmpchanlocs;' ...
                  'set(findobj(gcbf, ''tag'', ''ChanIndexEdit''), ''string'', vararg2str(ChanLabelIndex))']}, ...
        ... %5
        {'Style', 'text', 'string', 'Event type(s) to display and edit:'}, ...
        {'Style', 'edit', 'string', eventlist, 'tag', 'PatIDEventTypeEdit'}, ...
        {'Style', 'pushbutton', 'string', '...', ... 
                  'callback', ['tmpevent = EEG.event; [EventTypeIndex,EventTypeStr,EventTypeCell]=pop_chansel(unique({tmpevent.type})); clear tmpevent;' ...
                  'set(findobj(gcbf, ''tag'', ''PatIDEventTypeEdit''), ''string'', vararg2str(EventTypeCell))']}, ...
        ... %6
        {}, ...
        ... %7
        {'Style', 'text', 'string', 'Optional input key/val pairs:.'}, ...
        ... %8
        {'Style', 'edit'}, ...
        ... %9
        {}, ...
     }, ...
     'pophelp(''pop_VisEd'');', 'Select visual editing parameters -- pop_VisEd()' ...
     );
 
     if isempty(results);return;end
     
     DataType=results{1};
     ChanIndex=results{2};
     EventType=results{3};
     Options=results{4};
end


% return command
% -------------------------
if isempty(Options);
    com=sprintf('EEG = pop_VisEd( %s, %d, %s, {%s});', inputname(1), DataType, vararg2str(ChanIndex), EventType);
else
    com=sprintf('EEG = pop_VisEd( %s, %d, %s, {%s}, %s);', inputname(1), DataType, vararg2str(ChanIndex), EventType, Options);
end

% call command
% ------------
if isempty(Options);
    exec=sprintf('EEG = VisEd( %s, %d, %s, {%s});', inputname(1), DataType, vararg2str(ChanIndex), EventType);
else
    exec=sprintf('EEG = VisEd( %s, %d, %s, {%s}, %s);', inputname(1), DataType, vararg2str(ChanIndex), EventType, Options);
end
eval(exec);
