function res = length(this);

    tmp = struct(this);
    %if any(cellfun(@length, { tmp.EEG }) > 1)
    %    error('EEG structure in object with more than 1 element')
    %end
    try
        res = length(tmp.EEG);
    catch
        res = length(tmp);
        return;
    end
    if res > 1, error('EEG structure in object with more than 1 element'); end
    res = length(tmp);
