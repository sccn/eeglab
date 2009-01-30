%% Copyright (C) 1995, 1996, 1997  Friedrich Leisch
%%
%% This file is part of Octave.
%%
%% Octave is free software; you can redistribute it and/or modify it
%% under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2, or (at your option)
%% any later version.
%%
%% Octave is distributed in the hope that it will be useful, but
%% WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%% General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with Octave; see the file COPYING.  If not, write to the Free
%% Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%% 02110-1301, USA.

%% -*- texinfo -*-
%% @deftypefn {Function File} {} hurst (@var{x})
%% Estimate the Hurst parameter of sample @var{x} via the rescaled range
%% statistic.  If @var{x} is a matrix, the parameter is estimated for
%% every single column.
%% @end deftypefn

%% Author: FL <Friedrich.Leisch@ci.tuwien.ac.at>
%% Description: Estimate the Hurst parameter
%% Modified by Alois Schloegl <a.schloegl@ieee.org> 
%%  adapted for the use with Matlab: 
%%  added support for data with missing samples (encoded with NaN)   

function H = hurst (x)

  if (nargin ~= 1)
      fprintf(1,'usage: hurst (x)');
  end

  if (isscalar (x))
    error ('hurst: x must not be a scalar')
  elseif (isvector (x))
    x = reshape (x, length (x), 1);
  end

  [xr, xc] = size (x);

  x0 = center(x,1);

  xr = sum(~isnan(x0)); 
  x0(isnan(x0))=0;

  s = sqrt(sum(x0.*x0,1)./max(0,xr-1));	% std
  w = cumsum (x0,1);
  RS = (max(w) - min(w)) ./ s;
  H = log (RS) ./ log (xr);

