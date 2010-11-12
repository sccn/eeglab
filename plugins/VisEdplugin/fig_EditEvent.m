function g = fig_EditEvent(g, Latency, EventType, EventIndex, Proc)

switch Proc; %if strcmp(Proc, 'New');
    
    case 'New'
        
    if ~isfield(g, 'newindex');
        g.newindex=g.nevents+1;
    else
        g.newindex=g.newindex+1;
    end
    
    % Create new event.
    if isempty(g.events);
        g.events(1).latency=Latency;
    else
        g.events(length(g.events)+1).latency=Latency;
    end
    
    g.events(length(g.events)).type=EventType;
    g.events(length(g.events)).urevent=length(g.events);
    g.events(length(g.events)).proc='new';
    g.events(length(g.events)).index=g.newindex;
    if length(g.datasize)==3;
        g.events(length(g.events)).epoch=ceil(Latency/g.datasize(3));
    end
    
    if ~isfield(g, 'eventupdate');
        updateindex=1;
    else
        updateindex=length(g.eventupdate)+1;
    end
    
    g.eventupdate(updateindex).latency=Latency;
    g.eventupdate(updateindex).type=EventType;
    g.eventupdate(updateindex).proc='new';
    g.eventupdate(updateindex).index=g.newindex;
    if length(g.datasize)==3;
        g.eventupdate(updateindex).epoch=ceil(Latency/g.datasize(3));
    end

        
%end


case 'Delete'; %if strcmp(Proc, 'Delete');
    
    % log event update field.
    if ~isfield(g, 'eventupdate');
        updateindex=1;
    else
        updateindex=length(g.eventupdate)+1;
    end
    
    g.eventupdate(updateindex).latency=[];
    g.eventupdate(updateindex).type=[];
    g.eventupdate(updateindex).proc='clear';
    g.eventupdate(updateindex).index=g.events(EventIndex).index;
    
    % Clear SelEvent.
    g.events(EventIndex)=[];

%end


    case 'Edit';% if strcmp(Proc, 'Edit');
    
    
    % Set default "edittime" if edit time has not already been set for this
    % figure.
    if ~isfield(g, 'edittime');
        g.edittime=[-500 500];
    end
       
    % Define index of channel to display.
    chanind=round((g.chans-g.elecoffset)-(g.tmppos(1,2)*g.dispchans));
    
    
    % EditFig parameter UI goes here.
    % pop up window
    % -------------
    results=inputgui( ...
        {[1] [4 3] [4 3]}, ...
        {...
        ...1
        {'Style', 'text', 'string', 'Select edit figure parameters.', 'FontWeight', 'bold'}, ...
        ...2
        {'Style', 'text', 'string', 'Time around event to display (eg. -500 500):'}, ...
        {'Style', 'edit', 'string', num2str(g.edittime)}, ...
        ...3
        {'Style', 'text', 'string', 'Index of channel to be displayed:'}, ...
        {'Style', 'edit', 'string', num2str(chanind)}, ...
        }, ...
        'pophelp(''pop_EventEdit'');', 'event edit -- pop_EventEdit()' ...
        ); ...
        if isempty(results);return;end

    g.edittime = str2num(results{1});
    chanind    = str2num(results{2});
    
    
    % Check if EditTime extends outside of available data.
    % - if possible adjust window to available data.
    editpnts=g.edittime*(g.srate/1000);
    
    if Latency<=editpnts(1)*-1;
        EditStartPnt=1;
        MrkLat=Latency;
    else
        EditStartPnt=Latency+editpnts(1);
        MrkLat=abs(editpnts(1));
    end
    EditEndPnt=EditStartPnt+(editpnts(2)-editpnts(1));

    
    % Create EditData (data to be displayed in EditPlot.
    EditData=[];
    EditData=g.data(chanind, EditStartPnt:EditEndPnt);
    
    
    % if requested create EditData2 (overlay waveform).
    if ~isempty(g.data2);
        EditData2=g.data2(chanind, EditStartPnt:EditEndPnt);
    end


    % Set EditPlot axes parameters and target event limits.
    EditMin=min(EditData);
    EditMax=max(EditData);
    MrkMin=EditMin-((EditMax-EditMin)*.5);
    MrkMax=EditMax+((EditMax-EditMin)*.5);
    
    
    %this will be removed.
    EditTime=[EditStartPnt*(1000/g.srate):1000/g.srate:EditEndPnt*(1000/g.srate)];
    
    
    %create EditPlot fgure and call ginput.
    EditFig=figure;
    plot(EditTime,EditData, 'k');
    hold on;
    if ~isempty(g.data2);
        plot(EditTime,EditData2, 'r');
    end
        
    plot([EditTime(MrkLat) EditTime(MrkLat)+.001], [MrkMin MrkMax], '-b');
    axis([EditTime(1) EditTime(length(EditTime)) MrkMin MrkMax]);
    hold off;
    x=ginput;
    close(EditFig);
%    tmpLatency=g.events(EventIndex).latency;
    g.events(EventIndex).latency=round(x(length(x(:,1)),1)*(g.srate/1000));

    
    %Update "g.eventupdate" field.
    if ~isfield(g, 'eventupdate');
        updateindex=1;
    else
        updateindex=length(g.eventupdate)+1;
    end
    
    g.eventupdate(updateindex).latency=g.events(EventIndex).latency;
    g.eventupdate(updateindex).type=EventType;
    g.eventupdate(updateindex).proc='edit';
    g.eventupdate(updateindex).index=g.events(EventIndex).index;

end


% Create new eventtypes parameters if necessary.
if isfield(g, 'eventtypes');
    if ~any(strcmp(EventType, g.eventtypes));
        eventtypesN=length(g.eventtypes)+1;
        g.eventtypes{eventtypesN} = EventType;
        g.eventtypecolors{eventtypesN} = 'k';
        g.eventtypestyle{eventtypesN} = '-';
        g.eventtypewidths(eventtypesN) = 1;
    end
else
    eventtypesN=1;
    g.eventtypes{eventtypesN} = EventType;
    g.eventtypecolors{eventtypesN} = 'k';
    g.eventtypestyle{eventtypesN} = '-';
    g.eventtypewidths(eventtypesN) = 1;
    g.plotevent='on';
end
% Clear remaining display parameters.
if isfield(g, 'eventcolors');
    fields={'eventcolors', 'eventstyle', 'eventwidths', 'eventlatencies', 'eventlatencyend'};
    g=rmfield(g,fields);
end

if isempty(g.events);
    g.eventcolors=[];
    g.eventstyle=[];
    g.eventwidths=[];
    g.eventlatencies=[];
    g.eventlatencyend=[];    
else
    for i=1:length(g.events);
        eventtypeindex=find(strcmp(g.eventtypes,g.events(i).type));
        g.eventcolors{i}=g.eventtypecolors{eventtypeindex};
        g.eventstyle{i}=g.eventtypestyle{eventtypeindex};
        g.eventwidths(i)=g.eventtypewidths(eventtypeindex);
        g.eventlatencies(i)=g.events(i).latency;
        g.eventlatencyend(i)=g.events(i).latency+g.eventwidths(i);    
    end
end

g = rmfield(g, 'eventedit');

set(gcf, 'UserData', g);
eegplot('drawp', 0);

