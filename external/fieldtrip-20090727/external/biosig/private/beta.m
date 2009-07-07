%% Copyright (C) 2008 Alois Schloegl
%% $Id: beta.m,v 1.1 2009-07-07 02:23:48 arno Exp $
%% This function is part of BioSig http://biosig.sf.net 
%% Originally, it was part of Octave. It was modified for the use with FreeMat
%%
%% This is free software; you can redistribute it and/or modify it
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
%% @deftypefn {Mapping Function} {} beta (@var{a}, @var{b})
%% Return the Beta function,
%% @iftex
%% @tex
%% $$
%%  B (a, b) = {\Gamma (a) \Gamma (b) \over \Gamma (a + b)}.
%% $$
%% @end tex
%% @end iftex
%% @ifinfo
%%
%% @example
%% beta (a, b) = gamma (a) * gamma (b) / gamma (a + b).
%% @end example
%% @end ifinfo
%% @end deftypefn

%% Author: KH <Kurt.Hornik@wu-wien.ac.at>
%% Created: 13 June 1993
%% Adapted-By: jwe
%% Adapted for FreeMat by : AS <a.schloegl@ieee.org> 

function retval = beta (a, b)

  if (nargin ~= 2)
    usage ('beta (a, b)');
  end

  retval = exp (gammaln (a) + gammaln (b) - gammaln (a+b));

end
