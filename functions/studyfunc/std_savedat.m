% save computed data

function std_savedat( tmpfile, structure)

    delims = find( tmpfile == '.');
    structure.datafile = [ tmpfile(1:delims(end)-1) '.set' ];
    v = version;
    if v(1) < 7
        save('-mat', tmpfile, '-struct', 'structure');
    else
        save('-v6' , tmpfile, '-struct', 'structure');
    end;
