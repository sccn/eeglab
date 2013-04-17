% unique_bc - unique backward compatible with Matlab versions prior to 2013a

function [C,IA,IB] = unique_bc(A,varargin);

errorFlag = error_bc;

v = version;
indp = find(v == '.');
v = str2num(v(1:indp(2)-1));

if nargin > 2
    ind = strmatch('legacy', varargin);
    if ~isempty(ind)
        varargin(ind) = [];
    end;
end;

if v >= 7.14
    [C,IA,IB] = unique(A,varargin{:},'legacy');
    if warningFlag
        [C2,IA2,IB2] = unique(A,varargin{:});
        if warningFlag && (~isequal(C, C2) || ~isequal(IA, IA2) || ~isequal(IB, IB2))
            warning('backward compatibility issue with call to unique function');
        end;
    end;
else
    [C,IA,IB] = unique(A,varargin{:});
end;