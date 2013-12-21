function checkfunctionmatlab(func, toolbox)

if ismatlab
    if license('test',toolbox)
        p1 = fileparts(which('firls'));
        rmpath(p1);
        p2 = fileparts(which('firls'));
        if ~isempty(p2)
            error( [ 'Octave functions should not run on Matlab' 10 'Removing path to ' p1 '. Run your command again' ]);
        else
            addpath(p1);
            warning([ toolbox ' toolbox is absent or not in the path, using replacement functions' ]);
        end;
    else
        warning([ toolbox ' toolbox is absent or not in the path, using replacement functions' ]);
    end;
end;
