%% Copyright (C) 2000-2001 Paul Kienzle
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

%% -*- texinfo -*-
%% @deftypefn {Function File} {} unique (@var{x})
%%
%% Return the unique elements of @var{x}, sorted in ascending order.
%% If @var{x} is a row vector, return a row vector, but if @var{x}
%% is a column vector or a matrix return a column vector.
%%
%% @deftypefnx {Function File} {} unique (@var{A}, 'rows')
%%
%% Return the unique rows of @var{A}, sorted in ascending order.
%%
%% @deftypefnx {Function File} {[@var{y}, @var{i}, @var{j}] = } unique (@var{x})
%%
%% Return index vectors @var{i} and @var{j} such that @code{x(i)==y} and
%% @code{y(i)==x}.
%%
%% @end deftypefn
%% @seealso{union, intersect, setdiff, setxor, ismember}

%% modified by Alois Schloegl for the use with Matlab

function [y, i, j] = unique (x, r)

  if ( nargin < 1 || nargin > 2 || (nargin == 2 && ~strcmp(r,'rows')) )
    fprintf(1,'usage: unique (x) or unique (x, ''rows'')');
  end

  if (nargin == 1)
    n = prod(size(x));
  else
    n = size(x,1);
  end

  y = x; 
  if (n < 1)
    i = [];
    j = [];
    return
  elseif (n < 2)
    i = 1;
    j = 1;
    return
  end

  if isstr(x), y = toascii(y); end

  if nargin == 2
    [y, i] = sortrows(y);
    if (iscell(y))
      match = cellfun ('size',y(1:n-1,:),1) == cellfun ('size',y(2:n,:),1);
      idx = find(match);
      match(idx) = all(char(y(idx)) == char(y(idx+1)),2);
      match = all (match');
    else
      match = all( [ y(1:n-1,:) == y(2:n,:) ]' );
    end
    idx = find (match);
    y (idx, :) = [];
  else
    if (size(y,1) ~= 1) y = y(:); end
    [y, i] = sort(y);
    if (iscell(y))
      match = cellfun ('length',y(1:n-1)) == cellfun ('length',y(2:n));
      idx = find(match);
      match(idx) = all(char(y(idx)) == char(y(idx+1)),2);
    else
      match = [ y(1:n-1) == y(2:n) ];
    end
    idx = find (match);
    y (idx) = [];
  end

  %% I don't know why anyone would need reverse indices, but it
  %% was an interesting challenge.  I welcome cleaner solutions.
  if (nargout >= 3)
    j = i;
    j (i) = cumsum ( prepad ( ~match, n, 1 ) );
  end
  i (idx) = [];

  if isstr(x), y = char(y); end

end

%!assert(unique([1 1 2; 1 2 1; 1 1 2]),[1;2])
%!assert(unique([1 1 2; 1 0 1; 1 1 2],'rows'),[1 0 1; 1 1 2])
%!assert(unique([]),[])
%!assert(unique([1]),[1])
%!assert(unique([1 2]),[1 2])
%!assert(unique([1;2]),[1;2])
%!assert(unique([1,NaN,Inf,NaN,Inf]),[1,Inf,NaN,NaN])
%!assert(unique({'Foo','Bar','Foo'}),{'Bar','Foo'})
%!assert(unique({'Foo','Bar','FooBar'}),{'Bar','Foo','FooBar'})
