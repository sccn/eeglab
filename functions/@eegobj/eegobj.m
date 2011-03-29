function this = eegobj(EEG);

    if nargin == 1
        if isa(EEG, 'eegobj')
            this = EEG;
            return;
        end;
        for index = 1:length(EEG)
            TMP(index).EEG = EEG(index);
        end;
    else
        TMP.EEG = eeg_emptyset;
    end;
    this = class(TMP, 'eegobj');
