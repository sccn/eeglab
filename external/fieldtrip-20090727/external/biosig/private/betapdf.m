%% Copyright (C) 1995, 1996, 1997  Kurt Hornik
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
%% @deftypefn {Function File} {} betapdf (@var{x}, @var{a}, @var{b})
%% For each element of @var{x}, returns the PDF at @var{x} of the beta
%% distribution with parameters @var{a} and @var{b}.
%% @end deftypefn

%% Author: KH <Kurt.Hornik@wu-wien.ac.at>
%% Description: PDF of the Beta distribution

function pdf = betapdf (x, a, b)

  if (nargin ~= 3)
    error('usage: betapdf(x,a,b)');
  end
  
  if (~isscalar (a) || ~isscalar(b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
      error ('betapdf: x, a and b must be of common size or scalar');
    end
  end

  sz = size (x);
  pdf = zeros (sz);

  k = find (~(a > 0) | ~(b > 0) | isnan (x));
  if (any (k))
    pdf (k) = NaN;
  end

  k = find ((x > 0) & (x < 1) & (a > 0) & (b > 0));
  if (any (k))
    if (isscalar(a) && isscalar(b))
      pdf(k) = exp ((a - 1) .* log (x(k)) + (b - 1) .* log (1 - x(k))) ./ beta (a, b);
    else
      pdf(k) = exp ((a(k) - 1) .* log (x(k)) + (b(k) - 1) .* log (1 - x(k))) ./ beta (a(k), b(k));
    end
  end

end
