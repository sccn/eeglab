% ismember_bc - ismember backward compatible with Matlab versions prior to 2013a

function [C,IA] = ismember_bc(A,B,varargin);

errorFlag = error_bc;

v = version;
indp = find(v == '.');
v = str2num(v(1:indp(2)-1));
if v > 7.19, v = floor(v) + rem(v,1)/10; end;

if nargin > 2
    ind = strmatch('legacy', varargin);
    if ~isempty(ind)
        varargin(ind) = [];
    end;
end;

if v >= 7.14
    [C,IA] = ismember(A,B,varargin{:},'legacy');
    if errorFlag
        [C2,IA2] = ismember(A,B,varargin{:});
        if (~isequal(C, C2) || ~isequal(IA, IA2))
            warning('backward compatibility issue with call to ismember function');
        end;
    end;
else
    [C,IA] = ismember(A,B,varargin{:});
end;