function EEG = VisEd_UpdateEvents(EEG, g)

if isfield(g, 'eventupdate');
    for i=1:length(g.eventupdate);
        if strcmp(g.eventupdate(i).proc, 'new');
            eventindex=length(EEG.event)+1;
            EEG.event(eventindex).latency=g.eventupdate(i).latency;
            EEG.event(eventindex).type=g.eventupdate(i).type;
            if ndims(EEG.data)==3;
                EEG.event(eventindex).epoch=g.eventupdate(i).epoch;
            end
        end
        if strcmp(g.eventupdate(i).proc, 'clear');
            eventindex=g.eventupdate(i).index;
            EEG.event(eventindex).action='clear';
        end
        if strcmp(g.eventupdate(i).proc, 'edit');
            eventindex=g.eventupdate(i).index;
            EEG.event(eventindex).action='edit';
            EEG.event(eventindex).actlat=g.eventupdate(i).latency;
        end
    end
end

j=0;
for i=1:length(EEG.event);
    if isfield(EEG.event(i),'action');
        if strcmp(EEG.event(i).action,'edit');
            EEG.event(i).latency=EEG.event(i).actlat;
        end
        if strcmp(EEG.event(i).action,'clear');
            j=j+1;
            clearInd(j)=i;
        end
    end
end

if isfield(EEG.event, 'action');
    EEG.event=rmfield(EEG.event,'action');
end
if isfield(EEG.event, 'actlat');
    EEG.event=rmfield(EEG.event,'actlat');
end

if exist('clearInd');
    EEG.event(clearInd)=[];
end

%sort events.
if ~isempty(EEG.event);
    tmpevent  = EEG.event;
    eventorder=[1:length(EEG.event);[tmpevent.latency]]';
    eventorder=sortrows(eventorder,2);
    for i=1:length(EEG.event);
        TMP.event(i)=EEG.event(eventorder(i,1));
    end
else
    TMP.event=[];
end

rmfield(EEG,'event');
EEG.event=TMP.event;

EEG=eeg_checkset(EEG, 'eventconsistency');


%remove data segments if required.
if ~isempty(g.winrej);
    if EEG.trials==1;
        [EEG LASTCOM] = eeg_eegrej(EEG,eegplot2event(g.winrej, -1));
    else
        EEG.reject.rejmanual=zeros(1,EEG.trials);
        for i=1:length(g.winrej(:,1));
            EEG.reject.rejmanual((g.winrej(i,1)/EEG.pnts)+1)=1;
        end
        EEG = pop_rejepoch( EEG, EEG.reject.rejmanual, 0);
    end
end
        
if isfield(g,'eloc_file');
    if length(g.eloc_file(1).labels)>=4;
        if strmatch(g.eloc_file(1).labels(1:4),'comp');
            datoric=2;
        end
    end
end
if ~exist('datoric', 'var');
    datoric=1;
end

if isfield(g.eloc_file, 'badchan');
    switch datoric
        case 1
            for i=1:length(g.eloc_file);
                EEG.chanlocs(g.eloc_file(i).index).badchan=g.eloc_file(i).badchan;
            end
        case 2
            for i=1:length(g.eloc_file);
                EEG.reject.gcompreject(g.eloc_file(i).index)=g.eloc_file(i).badchan;
            end
    end
end


eeglab redraw

clear EEGTMP tmpcom;

