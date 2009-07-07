## Copyright (C) 2001 Paul Kienzle
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{n}, @var{d}]} = rat (@var{x}, @var{tol})
##
## Find a rational approximation to @var{x} within tolerance defined
## by @var{tol} using a continued fraction expansion. E.g,
##
## @example
##    rat(pi) = 3 + 1/(7 + 1/16) = 355/113
##    rat(e) = 3 + 1/(-4 + 1/(2 + 1/(5 + 1/(-2 + 1/(-7))))) = 1457/536
## @end example
##
## @end deftypefn
## @seealso{rats}

function [n,d] = rat(x,tol)

  if (nargin != [1,2] || nargout != 2)
    usage("[n,d] = rat(x,tol)");
  endif

  y = x(:);

  ## replace inf with 0 while calculating ratios
  y(isinf(y)) = 0;

  ## default norm
  if (nargin < 2)
    tol = 1e-6 * norm(y,1);
  endif

  ## First step in the approximation is the integer portion
  n = round(y);  # first element in the continued fraction
  d = ones(size(y));
  frac = y-n;
  lastn = ones(size(y));
  lastd = zeros(size(y));

  ## grab new factors until all continued fractions converge
  while (1)
    ## determine which fractions have not yet converged
    idx = find (abs(y-n./d) >= tol);
    if (isempty(idx)) break; endif

    ## grab the next step in the continued fraction
    flip = 1./frac(idx);
    step = round(flip); # next element in the continued fraction
    frac(idx) = flip-step;

    ## update the numerator/denominator
    nextn = n;
    nextd = d;
    n(idx) = n(idx).*step + lastn(idx);
    d(idx) = d(idx).*step + lastd(idx);
    lastn = nextn;
    lastd = nextd;
  endwhile

  ## move the minus sign to the top
  n = n.*sign(d);
  d = abs(d);

  ## return the same shape as you receive
  n = reshape(n, size(x));
  d = reshape(d, size(x));

  ## use 1/0 for Inf
  n(isinf(x)) = sign(x(isinf(x)));
  d(isinf(x)) = 0;

endfunction
