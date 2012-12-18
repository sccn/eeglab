function obj3 = bsxfun(f, obj1, obj2)

    if isa(obj2, 'mmo'), error('bsxfun does not work yet with the second argument being a memory mapped object'); end;
    data1 = memmapfile(obj1.dataFile, 'writable', obj1.writable, 'format', { 'single' obj1.dimensions 'x' });
    data2 = obj2;

    % make new memory mapped data file (blank)
    % --------------------------------
    newFileName = mmo.getnewfilename;
    fid = fopen(newFileName, 'w');
    s1  = size(obj1);
    s2  = size(obj2);
    ss.type = '()';
    ss.subs(1:length(s1)-1) = { ':' };
    if length(s1) == length(s2) && s1(end) == s2(end)
        error('bsxfun does not work on the last dimension');
    end;
    for index = 1:s1(end) % scan last dimension
        ss.subs{length(s1)} = index;
        tmpdata = bsxfun(f,subsref(data1.Data.x, ss), data2);
        fwrite(fid, tmpdata, 'float');
    end;
    fclose(fid);    
    
    % create object
    % -------------
    obj3 = mmo(newFileName, s1, true, obj1.transposed);
    obj3 = updateWorkspace(obj3);
