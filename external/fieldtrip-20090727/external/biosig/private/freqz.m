%% Copyright (C) 1996, 1997 John W. Eaton
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
%% Software Foundation, 59 Temple Place - Suite 330, Boston, MA
%% 02111-1307, USA.

%% -*- texinfo -*-
%% @deftypefn {Function File} {[@var{h}, @var{w}] =} freqz (@var{b}, @var{a}, @var{n}, "whole")
%% Return the complex frequency response @var{h} of the rational IIR filter
%% whose numerator and denominator coefficients are @var{b} and @var{a},
%% respectively.  The response is evaluated at @var{n} angular frequencies
%% between 0 and
%% @ifinfo
%%  2*pi.
%% @end ifinfo
%% @iftex
%% @tex
%%  $2\pi$.
%% @end tex
%% @end iftex
%%
%% @noindent
%% The output value @var{w} is a vector of the frequencies.
%%
%% If the fourth argument is omitted, the response is evaluated at
%% frequencies between 0 and
%% @ifinfo
%%  pi.
%% @end ifinfo
%% @iftex
%% @tex
%%  $\pi$.
%% @end tex
%% @end iftex
%%
%% If @var{n} is omitted, a value of 512 is assumed.
%%
%% If @var{a} is omitted, the denominator is assumed to be 1 (this
%% corresponds to a simple FIR filter).
%%
%% For fastest computation, @var{n} should factor into a small number of
%% small primes.
%%
%% @deftypefnx {Function File} {@var{h} =} freqz (@var{b}, @var{a}, @var{w})
%% Evaluate the response at the specific frequencies in the vector @var{w}.
%% The values for @var{w} are measured in radians.
%%
%% @deftypefnx {Function File} {[@dots{}] =} freqz (@dots{}, @var{Fs})
%% Return frequencies in Hz instead of radians assuming a sampling rate
%% @var{Fs}.  If you are evaluating the response at specific frequencies 
%% @var{w}, those frequencies should be requested in Hz rather than radians.
%%
%% @deftypefnx {Function File} {} freqz (@dots{})
%% Plot the pass band, stop band and phase response of @var{h} rather
%% than returning them.
%% @end deftypefn

%% Author: jwe ???

function [h_r, w_r] = freqz (b, a, n, region, Fs)

  if (nargin < 1 || nargin > 5)
    usage ("[h, w] = freqz (b, a, n [, \"whole\"] [, Fs])");
  elseif (nargin == 1)
    %% Response of an FIR filter.
    a = n = region = Fs = [];
  elseif (nargin == 2)
    %% Response of an IIR filter
    n = region = Fs = [];
  elseif (nargin == 3)
    region = Fs = [];
  elseif (nargin == 4)
    Fs = [];
    if (! isstr (region) && ! isempty (region))
      Fs = region; 
      region = [];
    endif
  endif

  if (isempty (a)) 
    a = 1; 
  endif
  if (isempty (n))
    n = 512; 
  endif
  if (isempty (region))
    if (isreal (b) && isreal (a))
      region = "half";
    else
      region = "whole";
    endif
  endif
  if (isempty (Fs)) 
    if (nargout == 0) 
      Fs = 2; 
    else 
      Fs = 2*pi; 
    endif
  endif

  la = length (a);
  a = reshape (a, 1, la);
  lb = length (b);
  b = reshape (b, 1, lb);
  k = max ([la, lb]);

  if (! isscalar (n))
    if (nargin == 4) %% Fs was specified
      w = 2*pi*n/Fs;
    else
      w = n;
    endif
    n = length (n);
    extent = 0;
  elseif (strcmp (region, "whole"))
    w = 2 * pi * (0:n-1) / n;
    extent = n;
  else
    w = pi * (0:n-1) / n;
    extent = 2 * n;
  endif

  if (length (b) == 1)
    if (length (a) == 1)
      hb = b * ones (1, n);
    else
      hb = b;
    endif
  elseif (extent >= k) 
    hb = fft (postpad (b, extent));
    hb = hb(1:n);
  else
    hb = polyval (postpad (b, k), exp (j*w));
  endif
  if (length (a) == 1)
    ha = a;
  elseif (extent >= k)
    ha = fft (postpad (a, extent));
    ha = ha(1:n);
  else
    ha = polyval (postpad (a, k), exp (j*w));
  endif
  h = hb ./ ha;
  w = Fs * w / (2*pi);

  if (nargout != 0), # return values and don't plot
    h_r = h;
    w_r = w;
  else             # plot and don't return values
    freqz_plot (w, h);
  end

%% endfunction
