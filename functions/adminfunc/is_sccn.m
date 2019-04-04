function bool = is_sccn
% is_sccn() - returns 1 if computer is located at SCCN (Swartz Center
% for computational Neuroscience) and 0 otherwise

bool = 0;
domnane = ' ';
try
    eval([ 'if isunix, [tmp domname] = unix(''hostname'');' ...
        'end;' ...
        'bool = findstr(domname(1:end-1), ''ucsd.edu'');' ...
        'if isempty(bool), bool=0; else bool=1; end;' ], '');
catch
end

