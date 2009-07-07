function keyvalcheck(arglist, varargin)

% KEYVALCHECK is a helper function for parsing optional key-value input pairs.
%
% Use as
%   keyvalcheck(argin, 'required',  {'key1', 'key2', ...})
%   keyvalcheck(argin, 'forbidden', {'key1', 'key2', ...})
%   keyvalcheck(argin, 'optional',  {'key1', 'key2', ...})
%
% See also KEYVAL

% Copyright (C) 2009, Robert Oostenveld
%
% $Log: not supported by cvs2svn $
% Revision 1.2  2009/06/15 14:23:26  roboos
% only check an option if specified (i.e. non-empty)
%
% Revision 1.1  2009/04/14 19:37:44  roboos
% new helper function, used in plot_xxx functions
%

required  = keyval('required', varargin);
forbidden = keyval('forbidden', varargin);
optional  = keyval('optional', varargin);

keys = arglist(1:2:end);
vals = arglist(2:2:end);

if numel(keys)~=numel(vals)
  error('optional input arguments should come in key-value pairs, i.e. there should be an even number');
end

if ~isempty(required)
  % only check if specified
  set = intersect(keys, required);
  if numel(set)~=numel(required)
    error('the required input argument ''%s'' was not specified', set{:});
  end
end

if ~isempty(forbidden)
  % only check if specified
  set = intersect(keys, forbidden);
  if numel(set)~=0
    error('the input argument ''%s'' is forbidden', set{:});
  end
end

if ~isempty(optional)
  % only check if specified
  set = setdiff(keys, optional);
  if numel(set)>0
    error('the input argument ''%s'' is forbidden', set{:});
  end
end
