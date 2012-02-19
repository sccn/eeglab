function EEG = VisEd(EEG, DataType, ChanIndex, EventType,varargin)

%COLLECT AND SET g STRUDCTURE FROM VARARGIN KEY/VAL PAIRS...
g = struct(varargin{:});

try g.quick_evtmk; catch, g.quick_evtmk = ''; end; 
try g.quick_evtrm; catch, g.quick_evtrm = 'off'; end; 

chans=eval(ChanIndex);

if ~isfield(EEG.chanlocs,'badchan')
    for i=1:EEG.nbchan;
        EEG.chanlocs(i).badchan=0;
    end
end
switch DataType

    case 1
        % Create data vector containing EEG channels selected in pop_VisEd
        data=EEG.data(chans,:,:);

        % Create VisEd.chan field containing channel information for channels in
        % data vector.
        
        for i=1:length(EEG.chanlocs);
            EEG.chanlocs(i).index=i;
        end
        VisEd.chan=EEG.chanlocs(chans);
    
    case 2
        eeglab_options; % changed from eeglaboptions 3/30/02 -sm
        if option_computeica
            data = EEG.icaact;
        else
            data = (EEG.icaweights*EEG.icasphere)*reshape(EEG.data, length(EEG.icaweights(1,:)), EEG.trials*EEG.pnts);
            data = reshape( data, size(data,1), EEG.pnts, EEG.trials);
        end
        
        tmpdata=data(chans,:,:);
        data=[];
        data=tmpdata;

        for i=1:length(chans);
            VisEd.chan(i).labels=sprintf('%s%s','comp',num2str(chans(i)));
            VisEd.chan(i).badchan=EEG.reject.gcompreject(chans(i));
            VisEd.chan(i).index=chans(i);
        end

    case 3
        %procedure for plotting EXG composite channels. eventually...
end

% Create VisEd.event field containing events selected in pop_VisEd.

j=0;
if isempty(EventType);
    VisEd.event = [];
else
    for i=1:length(EEG.event);
        if ~isempty(strmatch(EEG.event(i).type,EventType, 'exact'));
            event=EEG.event(i);
            event.index=i;
            event.proc='none';
            j=j+1;
            VisEd.event(j)=event;
        end
    end
end

icacomp=2;
nbpnts=EEG.pnts;

if EEG.trials > 1
    if icacomp == 1 macrorej  = 'EEG.reject.rejmanual';
        			macrorejE = 'EEG.reject.rejmanualE';
    else			macrorej  = 'EEG.reject.icarejmanual';
        			macrorejE = 'EEG.reject.icarejmanualE';
    end;
	colrej = EEG.reject.rejmanualcol;
	rej  = eval(macrorej);
	rejE = eval(macrorejE);

    if all(colrej == EEG.reject.rejmanualcol)
        oldrej = [];  % for manual rejection, old rejection are
        oldrejE = []; % the current rejection
    else
        oldrej  = eval(macrorej);
        oldrejE = eval(macrorejE);
    end;

    rejeegplottmp = trial2eegplot(  oldrej, oldrejE, nbpnts, min(colrej+0.15, [1 1 1]));
    if ~isempty(rejeegplottmp)
        rejeegplot = [ rejeegplottmp ];
    else
        rejeegplot = [];
    end
    rejeegplottmp = trial2eegplot(  rej, rejE, nbpnts, colrej);
    if ~isempty(rejeegplottmp)
        rejeegplot = [ rejeegplot; rejeegplottmp ];
    end
    
else
    rejeegplot=EEG.reject.rejmanual;
end


command=sprintf('%s','EEG=VisEd_UpdateEvents(EEG, g); EEG.saved = ''no'';');

% Call eegplot with variables set from pop-VisEd and ctrlselectcommand
% option set to { 'VisEd_ctrldowncom;' 'eegplot(''defmotioncom'', gcbf);' '' }.
eegplot(data, ...
              'eloc_file', VisEd.chan, ...
              'events', VisEd.event, ...
              'srate', EEG.srate, ...
              'winrej', rejeegplot, ...
              'butlabel', 'Update EEG structure', ...
              'command', command, ...
              'ctrlselectcommand',{ ['VisEd_ctrldowncom(EEG,''', g.quick_evtmk, ''',''', g.quick_evtrm ''');'] 'eegplot(''defmotioncom'', gcbf);' '' });
              % 'extselectcommand',{ 'VisEd_extdowncom;' 'eegplot(''defmotioncom'', gcbf);' '' } ...
        


    
