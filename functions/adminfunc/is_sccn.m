% is_sccn() - returns 1 if computer is located at SCCN (Swartz Center
% for computational Neuroscience) and 0 otherwise

function bool = is_sccn;
    
    bool = 0;
    domnane = ' ';
    try 
        eval([ 'if isunix, [tmp domname] = unix(''hostname -d'');' ...
               'end;' ...
               'bool = strcmpi(domname(1:end-1), ''ucsd.edu'');' ], '');
    catch,
    end;
    
