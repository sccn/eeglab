function res = ctranspose(obj,useconj);
    if nargin == 1
        useconj = 1;
    end;

    if length(obj.dimensions) > 2
        error('Cannot transpose array');
    end;
    
    % make new memory mapped data file
    % --------------------------------
    newFileName = mmo.getnewfilename;
    copyfile(obj.dataFile, newFileName);
    
    res = obj;
    res.dimensions = [ obj.dimensions(2) obj.dimensions(1) ];
    res.dataFile = newFileName;
    tmpMMO1 = memmapfile(obj.dataFile, 'writable', obj.writable, 'format', { 'single' obj.dimensions 'x' });
    tmpMMO2 = memmapfile(res.dataFile, 'writable',         true, 'format', { 'single' res.dimensions 'x' });
    
    % copy the data
    % -------------
    if length(obj.dimensions) == 1 || obj.dimensions(1) > obj.dimensions(2)
        for index = 1:size(obj,2)
            s.type = '()';
            s.subs = { ':' index };
            if useconj, tmpMMO2.Data.x(index,:) = conj(subsref(tmpMMO1.Data.x,s));
            else        tmpMMO2.Data.x(:,index) =      subsref(tmpMMO1.Data.x,s);
            end;
        end;
    else
        for index = 1:size(obj,1)
            s.type = '()';
            s.subs = { index ':' };
            if useconj, tmpMMO2.Data.x(:,index) = conj(subsref(tmpMMO1.Data.x,s));
            else        tmpMMO2.Data.x(:,index) =      subsref(tmpMMO1.Data.x,s);
            end;
        end;
    end;
    