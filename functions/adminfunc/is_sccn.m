% is_sccn() - returns 1 if computer is located at SCCN (Swartz Center
% for computational Neuroscience and 0 otherwise

function bool = is_sccn;
    
    domnane = '';
    if isunix, [tmp domname] = unix('hostname -d');
    else return;
    end;
    bool = strcmpi(domname(1:end-1), 'sccn.ucsd.edu');
    return
    