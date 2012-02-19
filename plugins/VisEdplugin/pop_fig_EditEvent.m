% pop_fig_EditEvent() - Set parameters for editing events from the eegplot figure window.
%
% Usage: 
%   >>  g = pop_fig_EditEvent(data, g, Latency, EventType, EventIndex, Proc);
%
% Inputs:
%   data       - EEG channel data being displayed in eegplot figure window.
%   g          - eegplot UserData
%   Latency    - data point of button press.
%   EventType  - Label of selected or new event.
%   EventIndex - Index of Event to be edited (0 if creating new event).
%   Proc       - procedure to use on selected event (New, Edit, Delete).
%
% Outputs:
%   EEG  - output dataset
%
% If there are no events near the time point of the button press the user
% is given the option of either entering a new event into the data or
% toggling a bad channel status. If the "Edit events" check box is selected
% the string in the "Event tupe" edit box will be the "type" of the new
% event (note that while the "Event editing procedure" popup menu is
% present in this UI it is only populated by "New"). If the "Toggle bad
% channel status" check box is selected the bad channel status of the 
% channel identified by the label in the "Channel selection" popup menu 
% will alternate (this alternation affects the EEG.chanlocs.badchan field
% by setting it to 0 or 1). 
%
% If there are events close to the time point of the button press the user
% is given the option of selecting among existing events (within +/-20
% points of the button press) using the "Event close to press" popup menu,
% then perform a procedure listed in the "Event editing procedure" (New,
% Edit, Delete). Note: In this case if a new event name is entered into
% the "Event type" edit box the only procedure that makes sense is "New",
% but the other options are still available in the "Event editing procedure"
% popup box and will produce errors if used.

%
% See also:
%   EEGLAB, eegplot, VisEd 

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



function [g,com]=pop_fig_EditEvent(g, Latency, EventType, EventIndex, Proc)
% the command output is a hidden output that does not have to
% be described in the header
com = ''; % this initialization ensure that the function will return something
          % if the user press the cancel button            
          % display help if not enough arguments
% ------------------------------------

if nargin < 1
	help pop_sig_EditEvent;
	return;
end;	


if nargin < 5
    
    ProcCell={'New', 'Delete'};
    
    if ~isfield(g.eventedit, 'SelEventStruct');
        tmp.event(1).index=0;
        tmp.event(1).Dist=0;
        tmp.event(1).type='User';
        tmp.event(1).latency=g.eventedit.PosLat;
        if length(g.datasize)==3;
            tmp.event(1).epoch=floor(tmp.event(1).latency/g.datasize(3));
        end
        
        if ~isempty(g.quick_evtmk);
            results = {1 g.quick_evtmk 1 0 {''}};
        else
        
            % pop up window
            % -------------
            
            results=inputgui( ...
                {[1] [1] [2 2] [2 2] [1] [1] [2 2] [1]}, ...
                {...
                ...1
                {'Style', 'text', 'string', 'Select editing parameters.', 'FontWeight', 'bold'}, ...
                ...2
                {'Style', 'checkbox', 'tag', 'EditEventCheck', 'string', 'Edit events:', 'value', 1, ...
                'callback', 'set(findobj(''tag'', ''MarkBadChanCheck''), ''value'', 0);' }, ...
                ...3
                {'Style', 'text', 'string', 'Event type.'}, ...
                {'Style', 'text', 'string', 'Event editing procedure.'}, ...
                ...4
                {'Style', 'edit', 'tag', 'SelEventTypeEdit', 'string', tmp.event(1).type}, ...
                {'Style', 'Popup', 'string', ProcCell{1}}, ...
                ...5
                {}, ...
                ...6
                {'Style', 'checkbox', 'tag', 'MarkBadChanCheck', 'string', 'Toggle bad channel status:', 'value', 0, ...
                'callback', 'set(findobj(''tag'', ''EditEventCheck''), ''value'', 0);' }, ...
                ...7
                {'Style', 'text', 'string', 'Channels selection:'}, ...
                {'Style', 'pushbutton', 'string', '...', 'tag', 'ChanLabelButton',...
                'callback', ['[ChanLabelIndex,ChanLabelStr,ChanLabelCell]=pop_chansel({g.eloc_file.labels});' ...
                'set(findobj(gcbf, ''tag'', ''ChanLabelEdit''), ''string'', vararg2str(ChanLabelCell))']}, ...
                ...8
                {'Style', 'edit', 'string', {g.eloc_file(g.eventedit.ChanIndex).labels} ,'tag', 'ChanLabelEdit'}, ...
                }, ...
                'pophelp(''pop_fig_EditEvent'');', 'event edit -- pop_fig_EditEvent()');%, [], 'return');
            %close;
            if isempty(results);return;end
    
        end
        
        if results{1}==1;
            Proc       = ProcCell{results{3}};
            
            Latency    = g.eventedit.PosLat;
            EventType  = results{2};
            EventIndex = 0;
        end
        if results{4}==1;
            if ~isfield(g.eloc_file, 'badchan');
                for i=1:length(g.eloc_file);
                    g.eloc_file(i).badchan=0;
                end
            end
            ChanLabelStr=results{5};
            if iscell(ChanLabelStr);
                g.eventedit.ChanLabelCell=ChanLabelStr;
            else
                g.eventedit.ChanLabelCell=eval(['{' ChanLabelStr '}']);
            end
            for i=1:length(g.eventedit.ChanLabelCell);
                g.eventedit.ChanIndex=strmatch(g.eventedit.ChanLabelCell{i},{g.eloc_file.labels},'exact');
                if g.eloc_file(g.eventedit.ChanIndex).badchan==0;
                    g.eloc_file(g.eventedit.ChanIndex).badchan=1;
                else
                    g.eloc_file(g.eventedit.ChanIndex).badchan=0;
                end
            end
            g = rmfield(g, 'eventedit');
            set(findobj('tag', 'EEGPLOT'), 'UserData', g);
            eegplot('drawp', 0);
            return
        end
    else
        
        for i=1:length(g.eventedit.SelEventStruct);
            tmpInd(i,1)=g.eventedit.SelEventStruct(i).dist;
            tmpInd(i,2)=i;
        end
        tmpSort=sortrows(tmpInd,1);
        clear tmpInd

        for i=1:length(tmpSort(:,1));
            tmp.event(i)=g.eventedit.SelEventStruct(tmpSort(i,2));
        end
        clear tmpSort

        if strcmp(g.quick_evtrm,'on');
            results = {'' 2 1};
            Proc       = ProcCell{results{2}};
            if strcmp(Proc,'New');
                Latency    = g.eventedit.PosLat;
                EventType  = results{1};
                EventIndex = 0;
            else
                Latency    = tmp.event(results{3}).latency;
                EventType  = tmp.event(results{3}).type;
                EventIndex = tmp.event(results{3}).index;
            end
        else

        
            % pop up window
            % -------------
            if nargin < 5
                
                results=inputgui( ...
                    {[1] [2 2] [2 2] [2 2] [2 2]}, ...
                    {...
                    ...1
                    {'Style', 'text', 'string', 'Select editing parameters.', 'FontWeight', 'bold'}, ...
                    ...2
                    {'Style', 'text', 'string', 'Event type.'}, ...
                    {'Style', 'text', 'string', 'Event editing procedure.'}, ...
                    ...3
                    {'Style', 'edit', 'tag', 'SelEventTypeEdit', 'string', tmp.event(1).type}, ...
                    {'Style', 'Popup', 'tag', 'EventProcPupup','string', ProcCell, 'Value', 2} ...
                    ...5
                    {'Style', 'text', 'string', 'Events close to press:'}, ...
                    {}, ...
                    ...6
                    {'Style', 'Popup', 'tag', 'SelEventTypePopup', 'string', {tmp.event.type}, ...
                    'callback', ['tmpeventtype=get(findobj(''tag'', ''SelEventTypePopup''),''string'');', ...
                    'cureventtype=tmpeventtype{get(findobj(''tag'', ''SelEventTypePopup''),''Value'')};', ...
                    'set(findobj(''tag'', ''SelEventTypeEdit''),''string'', cureventtype);']}, ...
                    {}, ...
                    }, ...
                    'pophelp(''pop_fig_EditEvent'');', 'event edit -- pop_fig_EditEvent()' ...
                    ); ...
                    if isempty(results);return;end
                
                Proc       = ProcCell{results{2}};
                
                if strcmp(Proc,'New');
                    Latency    = g.eventedit.PosLat;
                    EventType  = results{1};
                    EventIndex = 0;
                else
                    Latency    = tmp.event(results{3}).latency;
                    EventType  = tmp.event(results{3}).type;
                    EventIndex = tmp.event(results{3}).index;
                end
            end
        end
    end


end



% return the string command
% -------------------------
com = sprintf('g = pop_fig_EditEvent(%s, %s, %s, %s, %s);', inputname(1), vararg2str(Latency), vararg2str(EventType), vararg2str(EventIndex), vararg2str(Proc));

% call function "FFTStandard" on raw data.
% ---------------------------------------------------
g=fig_EditEvent(g, Latency, EventType, EventIndex, Proc);

return;
