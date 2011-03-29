function this = orderfields(this, vals);

    for index = 1:length(this)
        this(index).EEG = orderfields(this(index).EEG, vals);
    end;
