function b = iscellnumeric(C)
% Return 1 if all elements of cell array are numeric
%
% Tim Mullen, 2011, SCCN/INC, UCSD

b = all(cellfun(@(x) isnumeric(x),C));
