function varargout = uimenu(varargin);

if ismatlab
    varargout{1} = builtin('uimenu', varargin{:});
else
    varargout = { [] };
end;