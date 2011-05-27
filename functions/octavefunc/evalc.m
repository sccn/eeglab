## Copyright (C) 2009 Samuel W. Sirlin.  All rights reserved.
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the
## Free Software Foundation; either version 2, or (at your option) any
## later version.
##
## This is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
## for more details.

function r = evalc(x,y)
% function r = evalc(x,y)
% function r = evalc(x)
%
% this is just a cover function for eval, for compatibility
%
% 6/12/2005 sws

if nargin < 2
  r = eval(x);
else
  r = eval(x,y);
end
% end
