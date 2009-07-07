function dispsig(signalMatrix, range, titlestr);
%DISPSIG - deprecated!
%
% Please use icaplot instead.
%
%   See also ICAPLOT

% @(#)$Id: dispsig.m,v 1.1 2009-07-07 02:23:52 arno Exp $

fprintf('\nNote: DISPSIG is now deprecated! Please use ICAPLOT.\n');

if nargin < 3, titlestr = ''; end
if nargin < 2, range = 1:size(signalMatrix, 1); end

icaplot('dispsig',signalMatrix',0,range,range,titlestr);
