% support function for mmo object

function obj2 = unitaryopp(f, obj1, varargin);

    % make new memory mapped data file (blank)
    % ----------------------------------------
    newFileName = mmo.getnewfilename;
    data = memmapfile(obj1.dataFile, 'writable', obj1.writable, 'format', { 'single' obj1.dimensions 'x' });
    fid = fopen(newFileName, 'w');
    if obj1.dimensions(end) == 1, obj1.dimensions(end) = []; end
    ss.type = '()';
    ss.subs(1:length(obj1.dimensions)-1) = { ':' };
    for index = 1:obj1.dimensions(end)
        ss.subs{length(obj1.dimensions)} = index;
        tmpdata = f(subsref(data.Data.x, ss), varargin{:});
        fwrite(fid, tmpdata, 'float');
    end
    fclose(fid);    
    
    % copy the data
    % -------------
    obj2 = mmo(newFileName, obj1.dimensions, true, obj1.transposed);
    obj2 = updateWorkspace(obj2);
    
