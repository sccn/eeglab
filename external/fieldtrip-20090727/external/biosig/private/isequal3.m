%% Copyright (C) 2000, 2005, 2006, 2007 Paul Kienzle
%%
%% This file is part of Octave.
%%
%% Octave is free software; you can redistribute it and/or modify it
%% under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 3 of the License, or (at
%% your option) any later version.
%%
%% Octave is distributed in the hope that it will be useful, but
%% WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%% General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with Octave; see the file COPYING.  If not, see
%% <http://www.gnu.org/licenses/>.

%% Undocumented internal function.

%% -*- texinfo -*-
%% @deftypefn {Function File} {} isequal3 (@var{nans_compare_equal}, @var{x1}, @var{x2}, @dots{})
%% Return true if @var{x1}, @var{x2}, @dots{} are all equal and
%% @var{nans_compare_equal} evaluates to false.
%%
%% If @var{nans_compare_equal} evaluates to true, then assume NaN == NaN.
%% @seealso{isequal, isequalwithequalnans}
%% @end deftypefn

%% Modified by: William Poetra Yoga Hadisoeseno
%% modified by Alois Schloegl for the use with Matlab

%% Algorithm:
%%
%% 1. Determine the class of x
%% 2. If x is of the struct, cell, list or char class, for each
%%    argument after x, determine whether it has the same class
%%    and size as x.
%%    Otherwise, for each argument after x, verify that it is not
%%    of the struct, cell, list or char class, and that it has
%%    the same size as x.
%% 3. For each argument after x, compare it for equality with x:
%%    a. struct     compare each member by name, not by order (recursive)
%%    b. cell/list  compare each member by order (recursive)
%%    c. char       compare each member with strcmp
%%    d. <other>    compare each nonzero member, and assume NaN == NaN
%%                  if nans_compare_equal is nonzero.

function t = isequal3(nans_compare_equal, x, varargin)

  if (nargin < 2)
    print_usage ();
  end

  l_v = nargin - 2;

  %% Generic tests.

  %% All arguments must either be of the same class or they must be
  %% numeric values.
  t = (all (strcmp (class(x), cellfun (@class, varargin, 'UniformOutput', false))) || (isnumeric (x) && all (cellfun (@isnumeric, varargin, 'UniformOutput', true))));

  if (t)
    %% Test that everything has the same number of dimensions.
    s_x = size (x);
    s_v = cellfun (@size, varargin, 'UniformOutput', false);
    t = all (length (s_x) == cellfun (@length, s_v));
  end

  if (t)
    %% Test that everything is the same size since it has the same
    %% dimensionality.
    l_x = length (s_x);
    s_v = reshape ([s_v{:}], length (s_x), []);
    idx = 0;
    while (t && idx < l_x)
      idx=idx+1;
      t = all (s_x(idx) == s_v(idx, :));
    end
  end

  if (t)
    %% Check individual classes.
    if (isstruct (x))
      %% Test the number of fields.
      fn_x = fieldnames (x);
      l_fn_x = length (fn_x);
      fn_v = cellfun (@fieldnames, varargin, 'UniformOutput', false);
      t = all (l_fn_x == cellfun (@length, fn_v));

      %% Test that all the names are equal.
      idx = 0;
      s_fn_x = sort (fn_x);
      while (t && idx < l_v)
      	idx=idx+1;
	%% We'll allow the fieldnames to be in a different order.
	t = all (strcmp (s_fn_x, sort (fn_v{idx})));
      end

      idx = 0;
      while (t && idx < l_fn_x)
	%% Test that all field values are equal.
      	idx=idx+1;
	args = {nans_compare_equal, {x.(fn_x{idx})}};
	for argn = 1:l_v
	  args{argn+2} = {varargin{argn}.(fn_x{idx})};
	end
	%% Minimize function calls by calling for all the arguments at
	%% once.
        t = isequal3(args{:});
      end

    elseif (iscell (x))
      %% Check that each element of a cell is equal.
      l_x = numel (x);
      idx = 0;
      while (t && idx < l_x)
	idx=idx+1;
	args = {nans_compare_equal, x{idx}};
	for p = 1:l_v
	  args{p+2} = varargin{p}{idx};
	end
        t = isequal3 (args{:});
      end

    elseif (ischar (x))

      %% Sizes are equal already, so we can just make everything into a
      %% row and test the rows.
      for i = 1:l_v
	strings{i} = reshape (varargin{i}, 1, []);
      end
      t = all (strcmp (reshape (x, 1, []), strings));

    else
      %% Check the numeric types.

      if (issparse (x))
	f_x = spfind (x);
      else
	f_x = find (x);
      end
      l_f_x = length (f_x);
      x = x(f_x);
      for argn = 1:l_v
	y = varargin{argn};
	if (issparse (y))
          f_y = spfind (y);
	else
          f_y = find (y);
	end

	t = (l_f_x == length (f_y)) && all (f_x == f_y);
	if (~t)
          return;
	end

	y = y(f_y);
	m = (x == y);
	t = all (m);

	if (~t)
          if (nans_compare_equal)
            t = isnan (x(~m)) && isnan (y(~m));
          else
            return;
          end
	end
      end

    end
  end

end

%% test size and shape
%!assert(isequal3(0,[1,2,3,4],[1,2,3,4]), true)
%!assert(isequal3(0,[1;2;3;4],[1;2;3;4]), true)
%!assert(isequal3(0,[1,2,3,4],[1;2;3;4]), false)
%!assert(isequal3(0,[1,2,3,4],[1,2;3,4]), false)
%!assert(isequal3(0,[1,2,3,4],[1,3;2,4]), false)

%!test
%! A = 1:8;
%! B = reshape (A, 2, 2, 2);
%! assert (isequal3 (0, A, B), false);

%!test
%! A = reshape (1:8, 2, 2, 2);
%! B = A;
%! assert (isequal3 (0, A, B), true);

%!test
%! A = reshape (1:8, 2, 4);
%! B = reshape (A, 2, 2, 2);
%! assert (isequal3 (0, A, B), false);

%% test for equality
%!assert(isequal3(0,[1,2,3,4],[1,2,3,4]), true)
%!assert(isequal3(1,{1,2,NaN,4},{1,2,NaN,4}), true)
%!assert(isequal3(1,[1,2,NaN,4],[1,2,NaN,4]), true)
%!assert(isequal3(0,['a','b','c','d'],['a','b','c','d']), true)
%% Test multi-line strings
%!assert(isequal3(0,["test";"strings"],["test";"strings"],["test";"strings"]), true)
%% test for inequality
%!assert(isequal3(0,[1,2,3,4],[1;2;3;4]),false)
%!assert(isequal3(0,{1,2,3,4},[1,2,3,4]),false)
%!assert(isequal3(0,[1,2,3,4],{1,2,3,4}),false)
%!assert(isequal3(0,[1,2,NaN,4],[1,2,NaN,4]),false)
%!assert(isequal3(1,[1,2,NaN,4],[1,NaN,3,4]),false)
%!assert(isequal3(1,[1,2,NaN,4],[1,2,3,4]),false)
%!assert(isequal3(0,['a','b','c','d'],['a';'b';'c';'d']),false)
%!assert(isequal3(0,{'a','b','c','d'},{'a';'b';'c';'d'}),false)
%% test for equality (struct)
%!assert(isequal3(0,struct('a',1,'b',2),struct('a',1,'b',2)),true)
%!assert(isequal3(0,struct('a',1,'b',2),struct('a',1,'b',2),struct('a',1,'b',2)),true)
%!assert(isequal3(0,struct('a','abc','b',2),struct('a','abc','b',2)),true)
%!assert(isequal3(1,struct('a',NaN,'b',2),struct('a',NaN,'b',2),struct('a',NaN,'b',2)),true)
%% test for inequality (struct)
%!assert(isequal3(0,struct('a',NaN,'b',2),struct('a',NaN,'b',2),struct('a',NaN,'b',2)),false)

