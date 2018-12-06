function channelloc_to_tsv(EEG)

% From an EEG variable (i.e. EEG=pop_loadset(*.set), export the channel
% location as tsv file following the BIDS specification
%
% FORMAT channelloc_to_tvs(EEG)
%
% Author: Cyril Pernet - LIMO Team, University of Edinurgh

% list of labels for which we know it's not an EEG channel
known_labels = {'EXG','TRIG','ECG','EOG','VEOG','HEOG','EMG','MISC'};

% channel.tsv
for electrode = 1:size(EEG.chanlocs,2)
    ename{electrode}     = EEG.chanlocs(electrode).labels; 
    if contains(EEG.chanlocs(electrode).labels,known_labels)
        type{electrode} = EEG.chanlocs(electrode).labels;
        if contains(EEG.chanlocs(electrode).labels,'EOG')
            unit{electrode} = [num2str(char(181)) 'V'];
        elseif contains(EEG.chanlocs(electrode).labels,'ECG')
            unit{electrode} = 'mV';
        else
            unit{electrode} = ' ';
        end
    else
        type{electrode} = 'EEG';
        unit{electrode} = [num2str(char(181)) 'V']; % char(181) is mu in ASCII
    end
    sampling_frequency(electrode)  = EEG.srate;
    reference{electrode} = EEG.ref;
end

t = table(ename',unit',type',sampling_frequency',reference','VariableNames',{'name','units','type','sampling_reference','reference'});
channels_tsv_name = [EEG.filepath filesep EEG.filename(1:end-4) '_channels.tsv'];
writetable(t,channels_tsv_name,'FileType','text','Delimiter','\t');

