function events_to_tsv(EEG)

% From an EEG variable (i.e. EEG=pop_loadset(*.set), export the EEG.event
% as tsv file 
%
% FORMAT electrodes_to_tsv(EEG)
%
% Author: Cyril Pernet - LIMO Team, University of Edinurgh

for event=1:size(EEG.event,2)
    onset(event)        = EEG.event(event).latency;
    trial_type{event}   = EEG.event(event).urevent;
    if ~isfield(EEG.event(event),'duration')
        duration(event) = 0;
    else
        duration(event) = EEG.event(event).duration;
    end        
end

t = table(onset',duration',trial_type','VariableNames',{'onset','duration','trial_type'});
events_tsv_name = [EEG.filepath filesep EEG.filename(1:end-4) '_events.tsv'];
writetable(t,events_tsv_name,'FileType','text','Delimiter','\t');


