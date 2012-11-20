function res = permute(obj, dims);

if length(dims) > 3
    error('Max 3 dimensions for permutation');
elseif length(dims) == 2
    error('Permutation with 2 dimensions: use transpose instead');
end;

newFileName = mmo.getnewfilename;
copyfile(obj.dataFile, newFileName);

res = obj;
res.dimensions = obj.dimensions(dims);
res.dataFile = newFileName;
tmpMMO1 = memmapfile(obj.dataFile, 'writable', obj.writable, 'format', { 'single' obj.dimensions 'x' });
tmpMMO2 = memmapfile(res.dataFile, 'writable', true        , 'format', { 'single' res.dimensions 'x' });

% copy the data
% -------------
d = res.dimensions;
for i1 = 1:obj.dimensions(1)
    s.type = '()';
    s.subs = { i1 ':' ':' };
    tmpdata = squeeze(subsref(tmpMMO1.Data.x,s));

    if all(dims == [2 1 3]), tmpMMO2.Data.x( :,i1, :) = reshape(tmpdata , d(1),    1, d(3)); end;
    if all(dims == [3 1 2]), tmpMMO2.Data.x( :,i1, :) = reshape(tmpdata', d(1),    1, d(3)); end;

    if all(dims == [2 3 1]), tmpMMO2.Data.x( :, :,i1) = reshape(tmpdata , d(1), d(2),    1); end;
    if all(dims == [3 2 1]), tmpMMO2.Data.x( :, :,i1) = reshape(tmpdata', d(1), d(2),    1); end;

    if all(dims == [1 2 3]), tmpMMO2.Data.x(i1, :, :) = reshape(tmpdata ,    1, d(2), d(3)); end;
    if all(dims == [1 3 2]), tmpMMO2.Data.x(i1, :, :) = reshape(tmpdata',    1, d(2), d(3)); end;
end;

% slower versions below
% for i1 = 1:obj.dimensions(1)
%     for i2 = 1:obj.dimensions(3)
%         s.type = '()';
%         s.subs = { i1 ':' i2 };
%         tmpdata = subsref(tmpMMO1.Data.x,s);
%         if all(dims == [2 1 3]), tmpMMO2.Data.x( :,i1,i2) = tmpdata; end;
%         if all(dims == [2 3 1]), tmpMMO2.Data.x( :,i2,i1) = tmpdata; end;
%         
%         if all(dims == [1 2 3]), tmpMMO2.Data.x(i1, :,i2) = tmpdata; end;
%         if all(dims == [3 2 1]), tmpMMO2.Data.x(i2, :,i1) = tmpdata; end;
%         
%         if all(dims == [1 3 2]), tmpMMO2.Data.x(i1,i2, :) = tmpdata; end;
%         if all(dims == [3 1 2]), tmpMMO2.Data.x(i2,i1, :) = tmpdata; end;
%     end;
% end;
% 
% for i1 = 1:obj.dimensions(2)
%     for i2 = 1:obj.dimensions(3)
%         s.type = '()';
%         s.subs = { ':' i1 i2 };
%         tmpdata = subsref(tmpMMO1.Data.x,s);
%         if all(dims == [1 2 3]), tmpMMO2.Data.x( :,i1,i2) = tmpdata; end;
%         if all(dims == [1 3 2]), tmpMMO2.Data.x( :,i2,i1) = tmpdata; end;
%         if all(dims == [2 1 3]), tmpMMO2.Data.x(i1, :,i2) = tmpdata; end;
%         if all(dims == [2 3 1]), tmpMMO2.Data.x(i1,i2, :) = tmpdata; end;
%         if all(dims == [3 2 1]), tmpMMO2.Data.x(i2,i1, :) = tmpdata; end;
%         if all(dims == [3 1 2]), tmpMMO2.Data.x(i2, :,i1) = tmpdata; end;
%     end;
% end;
