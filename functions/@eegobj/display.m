function display(this);

    disp(inputname(1));
    if length(this) == 1
        struct(this.EEG)
    else
        TMP = struct(this);
        
        TMP2 = TMP(1).EEG;
        fieldorder = fieldnames(TMP2);
        for index = 2:length(TMP)
            TMP2(index) = orderfields(TMP(index).EEG, fieldorder);
        end;
        TMP2
    end;
    
