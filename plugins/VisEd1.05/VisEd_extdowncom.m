
ax1 = findobj('tag','backeeg','parent',gcbf);
tmppos = get(ax1, 'currentpoint');


% Store "UserData" and "temppos" variables to "g" structures.             
g=get(gcbf, 'UserData');
%g.tmppos=tmppos;
%g.nevents=length(EEG.event);

%if ~isfield(g.events, 'index');
%    for i=1:length(g.events);
%        g.events(i).index=i;
%        g.events(i).proc='none';
%    end
%end

% Define relative starting point for figure window latencies: "EEG.eventedit.WinStartPnt"
%if EEG.trials==1; % define WinStartPt. data point at wich current display window begins.
%    g.eventedit.WinStartPnt=g.time*EEG.srate;
%    g.eventedit.EpochIndex=1;
%else
%    g.eventedit.WindowIndex=g.time+1;
%    g.eventedit.WinStartPnt=g.eventedit.WindowIndex*EEG.pnts-(EEG.pnts-1);
%    g.eventedit.EpochIndex=ceil((g.tmppos(1,1)+g.eventedit.WinStartPnt)/EEG.pnts);
%end;

%g.eventedit.PosLat=round(g.tmppos(1,1)+g.eventedit.WinStartPnt);


% Identify selected channel.
% By default use the 'Eelec' tex display of eegplot.
tmpChanIndex=strmatch(get(findobj(gcf,'Tag','Eelec'),'string'),{g.eloc_file.labels},'exact');
if length(tmpChanIndex)==1;
    g.eventedit.ChanIndex=tmpChanIndex;
else
    % Otherwise calculate ChanIndex from tmppos.
    nWin=(g.chans-g.dispchans)+1;
    stepWin=1/nWin;
    if g.dispchans==g.chans;
        curWinStrt=0;
    else
        curWinStrt=floor((1-get(findobj('tag','eegslider'),'value'))/stepWin);
    end
    curWinEnd=curWinStrt+g.dispchans;

    YIndex=floor((tmppos(1,2)/(1/(g.dispchans+1)))-.5);
    g.eventedit.ChanIndex=(curWinEnd-YIndex);

    if g.eventedit.ChanIndex==0;
        g.eventedit.ChanIndex=1;
    end
    if g.eventedit.ChanIndex>EEG.nbchan;
        g.eventedit.ChanIndex=EEG.nbchan;
    end
end
clear tmpChanIndex

if ~isfield(g.eloc_file, 'badchan');
    for i=1:length(g.eloc_file);
        g.eloc_file(i).badchan=0;
    end
end

if g.eloc_file(g.eventedit.ChanIndex).badchan==0;
    g.eloc_file(g.eventedit.ChanIndex).badchan=1;
else
    g.eloc_file(g.eventedit.ChanIndex).badchan=0;
end

g = rmfield(g, 'eventedit');
set(findobj('tag', 'EEGPLOT'), 'UserData', g);
eegplot('drawp', 0);
return

% Check for event selection (if events already exist in dataset).
%if ~isempty(EEG.event);

    % Check for event selection (within +/-20 points of button press).
%    if isfield(g.eventedit, 'SelEventStruct');
%        g=rmfield(g.eventedit,'SelEventStruct');
%    end
%    j=0;
%    for i=1:length(g.events);
%        if abs(g.events(i).latency-g.eventedit.PosLat)<20;
%            j=j+1;
%           g.eventedit.SelEventStruct(j).index=i;
%            g.eventedit.SelEventStruct(j).dist=abs(g.events(i).latency-round(g.tmppos(1,1)+g.eventedit.WinStartPnt));
%            g.eventedit.SelEventStruct(j).type=g.events(i).type;
%           g.eventedit.SelEventStruct(j).latency=g.events(i).latency;
%        end
%    end
%end

% Call event edit UI.

%g = pop_fig_EditEvent(g);
