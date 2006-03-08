% save computed data

function std_savedat( tmpfile, structure)
    
    v = version;
    if v(1) < 7
        save('-mat', tmpfile, '-struct', 'structure');
    else
        save('-v6' , tmpfile, '-struct', 'structure');
    end;
