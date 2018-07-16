function obj3 = binaryopp(f, obj1, obj2)

    if isa(obj2, 'mmo'), tmpobj = obj2; obj2 = obj1; obj1 = tmpobj; clear tmpobj; end
    if ~isequal(size(obj1), size(obj2)) && prod(size(obj2)) ~= 1
        error('Matrix dimensions must agree.');
    end
    data1 = memmapfile(obj1.dataFile, 'writable', obj1.writable, 'format', { 'single' obj1.dimensions 'x' });
    if isa(obj2, 'mmo'), 
         data2 = memmapfile(obj2.dataFile, 'writable', obj2.writable, 'format', { 'single' obj2.dimensions 'x' });
    else data2 = obj2;
    end

    % make new memory mapped data file (blank)
    % --------------------------------
    newFileName = mmo.getnewfilename;
    fid = fopen(newFileName, 'w');
    s1  = size(obj1);
    ss.type = '()';
    ss.subs(1:length(s1)-1) = { ':' };
    for index = 1:s1(end)
        ss.subs{length(s1)} = index;
        if prod(size(data2)) == 1
             tmpdata = f(subsref(data1.Data.x, ss), data2);
        else tmpdata = f(subsref(data1.Data.x, ss), subsref(data2.Data.x, ss));
        end
        fwrite(fid, tmpdata, 'float');
    end
    fclose(fid);    
    
    % create object
    % -------------
    obj3 = mmo(newFileName, s1, true, obj1.transposed);
    obj3 = updateWorkspace(obj3);
