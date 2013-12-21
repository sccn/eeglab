% Copyright (C) 2006 Peter V. Lanspeary
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; If not, see <http://www.gnu.org/licenses/>.
%
% USAGE:
%   [spectra,freq] = pwelch(x,window,overlap,Nfft,Fs,
%                           range,plot_type,detrend,sloppy)
%     Estimate power spectral density of data "x" by the Welch (1967)
%     periodogram/FFT method.  All arguments except "x" are optional.
%         The data is divided into segments.  If "window" is a vector, each
%     segment has the same length as "window" and is multiplied by "window"
%     before (optional) zero-padding and calculation of its periodogram. If
%     "window" is a scalar, each segment has a length of "window" and a
%     Hamming window is used.
%         The spectral density is the mean of the periodograms, scaled so that
%     area under the spectrum is the same as the mean square of the
%     data.  This equivalence is supposed to be exact, but in practice there
%     is a mismatch of up to 0.5% when comparing area under a periodogram
%     with the mean square of the data.
%
%  [spectra,freq] = pwelch(x,y,window,overlap,Nfft,Fs,
%                          range,plot_type,detrend,sloppy,results)
%     Two-channel spectrum analyser.  Estimate power spectral density, cross-
%     spectral density, transfer function and/or coherence functions of time-
%     series input data "x" and output data "y" by the Welch (1967)
%     periodogram/FFT method.
%       pwelch treats the second argument as "y" if there is a control-string
%     argument "cross", "trans", "coher" or "ypower"; "power" does not force
%     the 2nd argument to be treated as "y".  All other arguments are
%     optional.  All spectra are returned in matrix "spectra". 
%
%  [spectra,Pxx_ci,freq] = pwelch(x,window,overlap,Nfft,Fs,conf,
%                                 range,plot_type,detrend,sloppy)
%  [spectra,Pxx_ci,freq] = pwelch(x,y,window,overlap,Nfft,Fs,conf,
%                                 range,plot_type,detrend,sloppy,results)
%     Estimates confidence intervals for the spectral density.
%     See Hint (7) below for compatibility options.  Confidence level "conf"
%     is the 6th or 7th numeric argument.  If "results" control-string 
%     arguments are used, one of them must be "power" when the "conf"
%     argument is present; pwelch can estimate confidence intervals only for
%     the power spectrum of the "x" data.  It does not know how to estimate
%     confidence intervals of the cross-power spectrum, transfer function or
%     coherence; if you can suggest a good method, please send a bug report.
%
% ARGUMENTS
% All but the first argument are optional and may be empty, except that
% the "results" argument may require the second argument to be "y".
%
% x           % [non-empty vector] system-input time-series data
% y           % [non-empty vector] system-output time-series data
%
% window      % [real vector] of window-function values between 0 and 1; the
%             %       data segment has the same length as the window.
%             %       Default window shape is Hamming.
%             % [integer scalar] length of each data segment.  The default
%             %       value is window=sqrt(length(x)) rounded up to the
%             %       nearest integer power of 2; see 'sloppy' argument.
%
% overlap     % [real scalar] segment overlap expressed as a multiple of
%             %       window or segment length.   0 <= overlap < 1,
%             %       The default is overlap=0.5 .
%
% Nfft        % [integer scalar] Length of FFT.  The default is the length
%             %       of the "window" vector or has the same value as the
%             %       scalar "window" argument.  If Nfft is larger than the
%             %       segment length, "seg_len", the data segment is padded
%             %       with "Nfft-seg_len" zeros.  The default is no padding.
%             %       Nfft values smaller than the length of the data
%             %       segment (or window) are ignored silently.
%
% Fs          % [real scalar] sampling frequency (Hertz); default=1.0
%
% conf        % [real scalar] confidence level between 0 and 1.  Confidence
%             %       intervals of the spectral density are estimated from
%             %       scatter in the periodograms and are returned as Pxx_ci.
%             %       Pxx_ci(:,1) is the lower bound of the confidence
%             %       interval and Pxx_ci(:,2) is the upper bound.  If there
%             %       are three return values, or conf is an empty matrix,
%             %       confidence intervals are calculated for conf=0.95 .
%             %       If conf is zero or is not given, confidence intervals
%             %       are not calculated. Confidence intervals can be 
%             %       obtained only for the power spectral density of x;
%             %       nothing else.
%
% CONTROL-STRING ARGUMENTS -- each of these arguments is a character string.
%   Control-string arguments must be after the other arguments but can be in
%   any order.
%  
% range     % 'half',  'onesided' : frequency range of the spectrum is
%           %       zero up to but not including Fs/2.  Power from
%           %       negative frequencies is added to the positive side of
%           %       the spectrum, but not at zero or Nyquist (Fs/2)
%           %       frequencies.  This keeps power equal in time and
%           %       spectral domains.  See reference [2].
%           % 'whole', 'twosided' : frequency range of the spectrum is
%           %       -Fs/2 to Fs/2, with negative frequencies
%           %       stored in "wrap around" order after the positive
%           %       frequencies; e.g. frequencies for a 10-point 'twosided'
%           %       spectrum are 0 0.1 0.2 0.3 0.4 0.5 -0.4 -0.3 -0.2 -0.1
%           % 'shift', 'centerdc' : same as 'whole' but with the first half
%           %       of the spectrum swapped with second half to put the
%           %       zero-frequency value in the middle. (See "help
%           %       fftshift".
%           % If data (x and y) are real, the default range is 'half',
%           % otherwise default range is 'whole'.
%
% plot_type % 'plot', 'semilogx', 'semilogy', 'loglog', 'squared' or 'db':
%           % specifies the type of plot.  The default is 'plot', which
%           % means linear-linear axes. 'squared' is the same as 'plot'.
%           % 'dB' plots "10*log10(psd)".  This argument is ignored and a
%           % spectrum is not plotted if the caller requires a returned
%           % value.
%
% detrend   % 'no-strip', 'none' -- do NOT remove mean value from the data
%           % 'short', 'mean' -- remove the mean value of each segment from
%           %                    each segment of the data.
%           % 'linear',       -- remove linear trend from each segment of
%           %                    the data.
%           % 'long-mean'     -- remove the mean value from the data before
%           %              splitting it into segments.  This is the default.
%
%   sloppy  % 'sloppy': FFT length is rounded up to the nearest integer
%           %       power of 2 by zero padding.  FFT length is adjusted
%           %       after addition of padding by explicit Nfft argument.
%           %       The default is to use exactly the FFT and window/
%           %       segment lengths specified in argument list.
%
%   results % specifies what results to return (in the order specified
%           %   and as many as desired).
%           % 'power' calculate power spectral density of "x"
%           % 'cross' calculate cross spectral density of "x" and "y"
%           % 'trans' calculate transfer function of a system with
%           %         input "x" and output "y"
%           % 'coher' calculate coherence function of "x" and "y"
%           % 'ypower' calculate power spectral density of "y"
%           %  The default is 'power', with argument "y" omitted. 
%
% RETURNED VALUES:
%   If return values are not required by the caller, the results are
%     plotted and nothing is returned.
%
%   spectra % [real-or-complex matrix] columns of the matrix contain results
%           %        in the same order as specified by "results" arguments.
%           %        Each column contains one of the result vectors.
%
%   Pxx_ci  % [real matrix] estimate of confidence interval for power
%           %        spectral density of x.  First column is the lower
%           %        bound.  Second column is the upper bound.
%
%   freq    % [real column vector] frequency values 
%
% HINTS
% 1) EMPTY ARGS:
%    if you don't want to use an optional argument you can leave it empty
%    by writing its value as [].
% 2) FOR BEGINNERS:
%    The profusion of arguments may make pwelch difficult to use, and an
%    unskilled user can easily produce a meaningless result or can easily
%    mis-interpret the result.  With real data "x" and sampling frequency
%    "Fs", the easiest and best way for a beginner to use pwelch is
%    probably "pwelch(x,[],[],[],Fs)".  Use the "window" argument to
%    control the length of the spectrum vector.  For real data and integer
%    scalar M, "pwelch(x,2*M,[],[],Fs)" gives an M+1 point spectrum.
%    Run "demo pwelch" (octave only).
% 3) WINDOWING FUNCTIONS:
%    Without a window function, sharp spectral peaks can have strong
%    sidelobes because the FFT of a data in a segment is in effect convolved
%    with a rectangular window.  A window function which tapers off
%    (gradually) at the ends produces much weaker sidelobes in the FFT.
%    Hann (hanning), hamming, bartlett, blackman, flattopwin etc are
%    available as separate Matlab/sigproc or Octave functions.  The sidelobes
%    of the Hann window have a roll-off rate of 60dB/decade of frequency.
%    The first sidelobe of the Hamming window is suppressed and is about 12dB
%    lower than the first Hann sidelobe, but the roll-off rate is only
%    20dB/decade.  You can inspect the FFT of a Hann window by plotting
%    "abs(fft(postpad(hanning(256),4096,0)))".
%    The default window is Hamming.
% 4) ZERO PADDING:
%    Zero-padding reduces the frequency step in the
%    spectrum, and produces an artificially smoothed spectrum.  For example,
%    "Nfft=2*length(window)" gives twice as many frequency values, but
%    adjacent PSD (power spectral density) values are not independent;
%    adjacent PSD values are independent if "Nfft=length(window)", which is
%    the default value of Nfft.
% 5) REMOVING MEAN FROM SIGNAL: 
%    If the mean is not removed from the signal there is a large spectral
%    peak at zero frequency and the sidelobes of this peak are likely to
%    swamp the rest of the spectrum.  For this reason, the default behaviour
%    is to remove the mean.  However, the matlab pwelch does not do this.
% 6) WARNING ON CONFIDENCE INTERVALS
%    Confidence intervals are obtained by measuring the sample variance of
%    the periodograms and assuming that the periodograms have a Gaussian
%    probability distribution.  This assumption is not accurate.  If, for
%    example, the data (x) is Gaussian, the periodogram has a Rayleigh
%    distribution.  However, the confidence intervals may still be  useful.
%    
% 7) COMPATIBILITY WITH Matlab R11, R12, etc
%    When used without the second data (y) argument, arguments are compatible
%    with the pwelch of Matlab R12, R13, R14, 2006a and 2006b except that
%     1) overlap is expressed as a multiple of window length ---
%        effect of overlap scales with window length
%     2) default values of length(window), Nfft and Fs are more sensible, and
%     3) Goertzel algorithm is not available so Nfft cannot be an array of
%        frequencies as in Matlab 2006b.
%    Pwelch has four persistent Matlab-compatibility levels.  Calling pwelch
%    with an empty first argument sets the order of arguments and defaults
%    specified above in the USAGE and ARGUMENTS section of this documentation.
%          prev_compat=pwelch([]);
%          [Pxx,f]=pwelch(x,window,overlap,Nfft,Fs,conf,...);
%    Calling pwelch with a single string argument (as described below) gives
%    compatibility with Matlab R11 or R12, or the R14 spectrum.welch
%    defaults.  The returned value is the PREVIOUS compatibility string.
%
%    Matlab R11:  For compatibility with the Matlab R11 pwelch:
%          prev_compat=pwelch('R11-');
%          [Pxx,f]=pwelch(x,Nfft,Fs,window,overlap,conf,range,units);
%          % units of overlap are "number of samples"
%          % defaults: Nfft=min(length(x),256), Fs=2*pi, length(window)=Nfft,
%          %           window=Hanning, do not detrend,
%          % N.B.  "Sloppy" is not available.
%
%    Matlab R12:  For compatibility with Matlab R12 to 2006a pwelch:
%          prev_compat=pwelch('R12+');
%          [Pxx,f]=pwelch(x,window,overlap,nfft,Fs,...);
%          % units of overlap are "number of samples"
%          % defaults: length(window)==length(x)/8, window=Hamming,
%          %           Nfft=max(256,NextPow2), Fs=2*pi, do not detrend
%          % NextPow2 is the next power of 2 greater than or equal to the
%          % window length. "Sloppy", "conf" are not available.  Default
%          % window length gives very poor amplitude resolution.
%
%    To adopt defaults of the Matlab R14 "spectrum.welch" spectrum object
%    associated "psd" method.
%          prev_compat=pwelch('psd');
%          [Pxx,f] = pwelch(x,window,overlap,Nfft,Fs,conf,...);
%          % overlap is expressed as a percentage of window length,
%          % defaults: length(window)==64, Nfft=max(256,NextPow2), Fs=2*pi
%          %           do not detrend
%          % NextPow2 is the next power of 2 greater than or equal to the
%          % window length. "Sloppy" is not available.
%          % Default window length gives coarse frequency resolution.
%
%
% REFERENCES
%  [1] Peter D. Welch (June 1967): 
%   "The use of fast Fourier transform for the estimation of power spectra:
%   a method based on time averaging over short, modified periodograms."
%   IEEE Transactions on Audio Electroacoustics, Vol AU-15(6), pp 70-73
%
%  [2] William H. Press and Saul A. Teukolsky and William T. Vetterling and
%               Brian P. Flannery",
%   "Numerical recipes in C, The art of scientific computing", 2nd edition,
%      Cambridge University Press, 2002 --- Section 13.7.
%  [3] Paul Kienzle (1999-2001): "pwelch", http://octave.sourceforge.net/


function [varargout] = pwelch(x,varargin)

checkfunctionmatlab('pwelch', 'signal_toolbox')

  %
  % COMPATIBILITY LEVEL
  % Argument positions and defaults depend on compatibility level selected
  % by calling pwelch without arguments or with a single string argument.
  %   native:      compatib=1; prev_compat=pwelch(); prev_compat=pwelch([]);
  %   matlab R11:  compatib=2; prev_compat=pwelch('R11-');
  %   matlab R12:  compatib=3; prev_compat=pwelch('R12+');
  %   spectrum.welch defaults:  compatib=4; prev_compat=pwelch('psd');
  % In each case, the returned value is the PREVIOUS compatibility string.
  %
compat_str = {[]; 'R11-'; 'R12+'; 'psd'};
persistent compatib;
if ( isempty(compatib) || compatib<=0 || compatib>4 )
  % legal values are 1, 2, 3, 4
  compatib = 1;
end
if ( nargin <= 0 )
  error( 'pwelch: Need at least 1 arg. Use "help pwelch".' );
elseif ( nargin==1 && (ischar(x) || isempty(x)) )
  varargout{1} = compat_str{compatib};
  if ( isempty(x) ) % native
    compatib = 1;
  elseif ( strcmp(x,'R11-') )
    compatib = 2;
  elseif ( strcmp(x,'R12+') )
    compatib = 3;
  elseif ( strcmp(x,'psd') )
    compatib = 4;
  else
    error( 'pwelch: compatibility arg must be empty, R11-, R12+ or psd' );
  end
  % return
%
% Check fixed argument
elseif ( isempty(x) || ~isvector(x) )
  error( 'pwelch: arg 1 (x) must be vector.' );
else
  %  force x to be COLUMN vector
  if ( size(x,1)==1 )
    x=x(:);
  end
  %
  % Look through all args to check if  cross PSD, transfer function or
  % coherence is required.  If yes, the second arg is data vector "y".
  arg2_is_y = 0;
  x_len = length(x);
  nvarargin = length(varargin);
  for iarg=1:nvarargin
    arg = varargin{iarg};
    if ( ~isempty(arg) && ischar(arg) && ...
         ( strcmp(arg,'cross') || strcmp(arg,'trans') || ...
           strcmp(arg,'coher') || strcmp(arg,'ypower') ))
      % OK. Need "y". Grab it from 2nd arg.
      arg = varargin{1};
      if ( nargin<2 || isempty(arg) || ~isvector(arg) || length(arg)~=x_len )
        error( 'pwelch: arg 2 (y) must be vector, same length as x.' );
      end
      % force  COLUMN vector
      y = varargin{1}(:);
      arg2_is_y = 1;
      break;
    end
  end
  %
  % COMPATIBILITY
  % To select default argument values, "compatib" is used as an array index.
  % Index values are   1=native,  2=R11,  3=R12,  4=spectrum.welch
  %
  %  argument positions:
  %  arg_posn = varargin index of window, overlap, Nfft, Fs and conf
  %             args respectively, a value of zero ==>> arg does not exist
  arg_posn = [1 2 3 4 5;  % native
              3 4 1 2 5;  % Matlab R11- pwelch
              1 2 3 4 0;  % Matlab R12+ pwelch
              1 2 3 4 5]; % spectrum.welch defaults
  arg_posn  = arg_posn(compatib,:) + arg2_is_y;
  %
  %  SPECIFY SOME DEFAULT VALUES for (not all) optional arguments
  %  Use compatib as array index.
  %  Fs = sampling frequency
  Fs        = [ 1.0 2*pi 2*pi 2*pi ];
  Fs        = Fs(compatib);            
  %  plot_type: 1='plot'|'squared'; 5='db'|'dB'
  plot_type = [ 1 5 5 5 ];
  plot_type = plot_type(compatib);
  %  rm_mean: 3='long-mean'; 0='no-strip'|'none'
  rm_mean   = [ 3 0 0 0 ];
  rm_mean   = rm_mean(compatib);
  % use max_overlap=x_len-1 because seg_len is not available yet
  % units of overlap are different for each version:
  %    fraction, samples, or percent
  max_overlap = [ 0.95 x_len-1 x_len-1 95];
  max_overlap = max_overlap(compatib);
  % default confidence interval
  %  if there are more than 2 return values and if there is a "conf" arg
  conf      = 0.95 * (nargout>2) * (arg_posn(5)>0);
  %
  is_win    = 0;    % =0 means valid window arg is not provided yet
  Nfft      = [];   % default depends on segment length
  overlap   = [];   % WARNING: units can be #samples, fraction or percentage
  range     = ~isreal(x) || ( arg2_is_y && ~isreal(y) );
  is_sloppy = 0;
  n_results = 0;
  do_power  = 0;   
  do_cross  = 0;
  do_trans  = 0;
  do_coher  = 0;
  do_ypower = 0;
%
%  DECODE AND CHECK OPTIONAL ARGUMENTS
  end_numeric_args = 0;
  for iarg = 1+arg2_is_y:nvarargin
    arg = varargin{iarg};
    if ( ischar(arg) )
      % first string arg ==> no more numeric args
      % non-string args cannot follow a string arg
      end_numeric_args = 1;
      %
      % decode control-string arguments
      if ( strcmp(arg,'sloppy') )
        is_sloppy = ~is_win || is_win==1;
      elseif ( strcmp(arg,'plot') || strcmp(arg,'squared') )
        plot_type = 1;
      elseif ( strcmp(arg,'semilogx') )
        plot_type = 2;
      elseif ( strcmp(arg,'semilogy') )
        plot_type = 3;
      elseif ( strcmp(arg,'loglog') )
        plot_type = 4;
      elseif ( strcmp(arg,'db') || strcmp(arg,'dB') )
        plot_type = 5;
      elseif ( strcmp(arg,'half') || strcmp(arg,'onesided') )
        range = 0;
      elseif ( strcmp(arg,'whole') || strcmp(arg,'twosided') )
        range = 1;
      elseif ( strcmp(arg,'shift') || strcmp(arg,'centerdc') )
        range = 2;
      elseif ( strcmp(arg,'long-mean') )
        rm_mean = 3;
      elseif ( strcmp(arg,'linear') )
        rm_mean = 2;
      elseif ( strcmp(arg,'short') || strcmp(arg,'mean') )
        rm_mean = 1;
      elseif ( strcmp(arg,'no-strip') || strcmp(arg,'none') )
        rm_mean = 0;
      elseif ( strcmp(arg, 'power' ) )
        if ( ~do_power )
          n_results = n_results+1;
          do_power = n_results;
        end
      elseif ( strcmp(arg, 'cross' ) )
        if ( ~do_cross )
          n_results = n_results+1;
          do_cross = n_results;
        end
      elseif ( strcmp(arg, 'trans' ) )
        if ( ~do_trans )
          n_results = n_results+1;
          do_trans = n_results;
        end
      elseif ( strcmp(arg, 'coher' ) )
        if ( ~do_coher )
          n_results = n_results+1;
          do_coher = n_results;
        end
      elseif ( strcmp(arg, 'ypower' ) )
        if ( ~do_ypower )
          n_results = n_results+1;
          do_ypower = n_results;
        end
      else
        error( 'pwelch: string arg %d illegal value: %s', iarg+1, arg ); 
      end
      % end of processing string args
      %
    elseif ( end_numeric_args )
      if ( ~isempty(arg) )
        % found non-string arg after a string arg ... oops
        error( 'pwelch: control arg must be string' );
      end
    %
    % first 4 optional arguments are numeric -- in fixed order
    %
    % deal with "Fs" and "conf" first because empty arg is a special default
    % -- "Fs" arg -- sampling frequency
    elseif ( iarg == arg_posn(4) )
      if ( isempty(arg) )
        Fs = 1;
      elseif ( ~isscalar(arg) || ~isreal(arg) || arg<0 )
        error( 'pwelch: arg %d (Fs) must be real scalar >0', iarg+1 );
      else
        Fs = arg;
      end
    %
    %  -- "conf" arg -- confidence level
    %    guard against the "it cannot happen" iarg==0
    elseif ( arg_posn(5) && iarg == arg_posn(5) )
      if ( isempty(arg) )
        conf = 0.95;
      elseif ( ~isscalar(arg) || ~isreal(arg) || arg < 0.0 || arg >= 1.0 )
        error( 'pwelch: arg %d (conf) must be real scalar, >=0, <1',iarg+1 );
      else
        conf = arg;
      end
    %
    % skip all empty args from this point onward
    elseif ( isempty(arg) )
      1;
    %
    %  -- "window" arg -- window function
    elseif ( iarg == arg_posn(1) )
      if ( isscalar(arg) )
        is_win = 1;
      elseif ( isvector(arg) )
        is_win = length(arg);
        if ( size(arg,2)>1 )  % vector must be COLUMN vector
          arg = arg(:);
        end
      else 
        is_win = 0;
      end
      if ( ~is_win )
        error( 'pwelch: arg %d (window) must be scalar or vector', iarg+1 );
      elseif ( is_win==1 && ( ~isreal(arg) || fix(arg)~=arg || arg<=3 ) )
        error( 'pwelch: arg %d (window) must be integer >3', iarg+1 );
      elseif ( is_win>1 && ( ~isreal(arg) || any(arg<0) ) )
        error( 'pwelch: arg %d (window) vector must be real and >=0',iarg+1);
      end
      window = arg;
      is_sloppy = 0;
    %
    % -- "overlap" arg -- segment overlap
    elseif ( iarg == arg_posn(2) )
      if (~isscalar(arg) || ~isreal(arg) || arg<0 || arg>max_overlap )
        error( 'pwelch: arg %d (overlap) must be real from 0 to %f', ...
               iarg+1, max_overlap );
      end
      overlap = arg;
    %
    % -- "Nfft" arg -- FFT length
    elseif ( iarg == arg_posn(3) )
      if ( ~isscalar(arg) || ~isreal(arg) || fix(arg)~=arg || arg<0 )
        error( 'pwelch: arg %d (Nfft) must be integer >=0', iarg+1 );
      end
      Nfft = arg;
    %
    else
      error( 'pwelch: arg %d  must be string', iarg+1 );
    end
  end
  if ( conf>0 && (n_results && ~do_power ) )
   error('pwelch: can give confidence interval for x power spectrum only' );
    end
  %
  % end DECODE AND CHECK OPTIONAL ARGUMENTS.
  %
  % SETUP REMAINING PARAMETERS
  % default action is to calculate power spectrum only
  if ( ~n_results )
    n_results = 1;
    do_power = 1;
  end
  need_Pxx = do_power || do_trans || do_coher;
  need_Pxy = do_cross || do_trans || do_coher;
  need_Pyy = do_coher || do_ypower;
  log_two = log(2);
  nearly_one = 0.99999999999;
  %
  % compatibility-options
  % provides exact compatibility with Matlab R11 or R12
  %
  % Matlab R11 compatibility
  if ( compatib==2 )
    if ( isempty(Nfft) )
      Nfft = min( 256, x_len );
    end
    if ( is_win > 1 )
      seg_len = min( length(window), Nfft );
      window = window(1:seg_len);
    else
      if ( is_win )
        % window arg is scalar
        seg_len = window;
      else
        seg_len = Nfft;
      end
      % make Hann window (don't depend on sigproc)
      xx = seg_len - 1;
      window = 0.5 - 0.5 * cos( (2*pi/xx)*[0:xx].' );
    end
  %
  % Matlab R12 compatibility
  elseif ( compatib==3 )
    if ( is_win > 1 )
      % window arg provides window function
      seg_len = length(window);
    else
      % window arg does not provide window function; use Hamming
      if ( is_win )
        % window arg is scalar
        seg_len = window;
      else
        % window arg not available; use R12 default, 8 windows
        % ignore overlap arg; use overlap=50% -- only choice that makes sense
        % this is the magic formula for 8 segments with 50% overlap
        seg_len = fix( (x_len-3)*2/9 );
      end
      % make Hamming window (don't depend on sigproc)
      xx = seg_len - 1;
      window = 0.54 - 0.46 * cos( (2*pi/xx)*[0:xx].' );
      end
    if ( isempty(Nfft) )
      Nfft = max( 256, 2^ceil(log(seg_len)*nearly_one/log_two) );
    end
  %
  % Matlab R14 psd(spectrum.welch) defaults
  elseif ( compatib==4 )
    if ( is_win > 1 )
      % window arg provides window function
      seg_len = length(window);
    else
      % window arg does not provide window function; use Hamming
      if ( is_win )
        % window arg is scalar
        seg_len = window;
      else
        % window arg not available; use default seg_len = 64
        seg_len = 64;
      end
      % make Hamming window (don't depend on sigproc)
      xx = seg_len - 1;
      window = 0.54 - 0.46 * cos( (2*pi/xx)*[0:xx].' );
      end
    % Now we know segment length,
    % so we can set default overlap as number of samples
    if ( ~isempty(overlap) )
      overlap = fix(seg_len * overlap / 100 );
      end
    if ( isempty(Nfft) )
      Nfft = max( 256, 2^ceil(log(seg_len)*nearly_one/log_two) );
    end
  %
  % default compatibility level
  else %if ( compatib==1 )
    % calculate/adjust segment length, window function
    if ( is_win > 1 )
      % window arg provides window function
      seg_len = length(window);
    else
      % window arg does not provide window function; use Hamming
      if ( is_win )       % window arg is scalar
        seg_len = window;
      else
        % window arg not available; use default length:
        % = sqrt(length(x)) rounded up to nearest integer power of 2
        if ( isempty(overlap) )
          overlap=0.5;
        end
        seg_len = 2 ^ ceil( log(sqrt(x_len/(1-overlap)))*nearly_one/log_two );
      end
      % make Hamming window (don't depend on sigproc)
      xx = seg_len - 1;
      window = 0.54 - 0.46 * cos( (2*pi/xx)*[0:xx].' );
      end
    % Now we know segment length,
    % so we can set default overlap as number of samples
    if ( ~isempty(overlap) )
      overlap = fix(seg_len * overlap);
      end
    %
    % calculate FFT length
    if ( isempty(Nfft) )
      Nfft = seg_len;
    end
    if ( is_sloppy )
      Nfft = 2 ^ ceil( log(Nfft) * nearly_one / log_two );
    end
  end
  % end of compatibility options
  %
  % minimum FFT length is seg_len
  Nfft = max( Nfft, seg_len );
  % Mean square of window is required for normalising PSD amplitude.
  win_meansq = (window.' * window) / seg_len;
  %
  % Set default or check overlap.
  if ( isempty(overlap) )
    overlap = fix(seg_len /2);
  elseif ( overlap >= seg_len )
    error( 'pwelch: arg (overlap=%d) too big. Must be <length(window)=%d',...
           overlap, seg_len );
  end
  %
  % Pad data with zeros if shorter than segment. This should not happen.
  if ( x_len < seg_len )
    x = [x; zeros(seg_len-x_len,1)];
    if ( arg2_is_y )
      y = [y; zeros(seg_len-x_len,1)];
      end
    x_len = seg_len;
    end
  % end SETUP REMAINING PARAMETERS
  %
  % 
  % MAIN CALCULATIONS
  % Remove mean from the data
  if ( rm_mean == 3 )
    n_ffts = max( 0, fix( (x_len-seg_len)/(seg_len-overlap) ) ) + 1;
    x_len  = min( x_len, (seg_len-overlap)*(n_ffts-1)+seg_len );
    if ( need_Pxx || need_Pxy )
      x = x - sum( x(1:x_len) ) / x_len;
    end
    if ( arg2_is_y || need_Pxy)
      y = y - sum( y(1:x_len) ) / x_len;
    end
  end
  %
  % Calculate and accumulate periodograms
  %   xx and yy are padded data segments
  %   Pxx, Pyy, Pyy are periodogram sums, Vxx is for confidence interval
  xx = zeros(Nfft,1);
  yy = xx;
  Pxx = xx;
  Pxy = xx;
  Pyy = xx;
  if ( conf>0 )
    Vxx = xx;
  else
    Vxx = [];
  end
  n_ffts = 0;
  for start_seg = [1:seg_len-overlap:x_len-seg_len+1]
    end_seg = start_seg+seg_len-1;
    % Don't truncate/remove the zero padding in xx and yy
    if ( need_Pxx || need_Pxy )
      if ( rm_mean==1 ) % remove mean from segment
        xx(1:seg_len) = window .* ( ...
          x(start_seg:end_seg) - sum(x(start_seg:end_seg)) / seg_len);
      elseif ( rm_mean == 2 ) % remove linear trend from segment
        xx(1:seg_len) = window .* detrend( x(start_seg:end_seg) );
      else % rm_mean==0 or 3
        xx(1:seg_len) = window .* x(start_seg:end_seg);
      end
      fft_x = fft(xx);
    end
    if ( need_Pxy || need_Pyy )
      if ( rm_mean==1 ) % remove mean from segment
        yy(1:seg_len) = window .* ( ...
          y(start_seg:end_seg) - sum(y(start_seg:end_seg)) / seg_len);
      elseif ( rm_mean == 2 ) % remove linear trend from segment
        yy(1:seg_len) = window .* detrend( y(start_seg:end_seg) );
      else % rm_mean==0 or 3
        yy(1:seg_len) = window .* y(start_seg:end_seg);
      end
      fft_y = fft(yy);
    end
    if ( need_Pxx )
      % force Pxx to be real; pgram = periodogram
      pgram = real(fft_x .* conj(fft_x));
      Pxx = Pxx + pgram;
      % sum of squared periodograms is required for confidence interval
      if ( conf>0 )
        Vxx = Vxx + pgram .^2;
        end
    end
    if ( need_Pxy )
      % Pxy (cross power spectrum) is complex. Do not force to be real.
      Pxy = Pxy + fft_y .* conj(fft_x);
    end
    if ( need_Pyy )
      % force Pyy to be real
      Pyy = Pyy + real(fft_y .* conj(fft_y));
    end
    n_ffts = n_ffts +1;
  end
  %
  % Calculate confidence interval
  %    -- incorrectly assumes that the periodogram has Gaussian probability
  %       distribution (actually, it has a single-sided (e.g. exponential)
  %       distribution.
  % Sample variance of periodograms is (Vxx-Pxx.^2/n_ffts)/(n_ffts-1).
  %    This method of calculating variance is more susceptible to round-off
  %  error, but is quicker, and for double-precision arithmetic and the
  %  inherently noisy periodogram (variance==mean^2), it should be OK.
  if ( conf>0 && need_Pxx )
    if ( n_ffts<2 )
      Vxx = zeros(Nfft,1);
    else
      % Should use student distribution here (for unknown variance), but tinv
      % is not a core Matlab function (is in statistics toolbox. Grrr)
      Vxx = (erfinv(conf)*sqrt(2*n_ffts/(n_ffts-1))) * sqrt(Vxx-Pxx.^2/n_ffts);
    end
  end
  %
  % Convert two-sided spectra to one-sided spectra (if range == 0).
  % For one-sided spectra, contributions from negative frequencies are added 
  % to the positive side of the spectrum -- but not at zero or Nyquist
  % (half sampling) frequencies.  This keeps power equal in time and spectral
  % domains, as required by Parseval theorem.
  %
  if ( range == 0 )
    if ( ~ rem(Nfft,2) )    % one-sided, Nfft is even
      psd_len = Nfft/2+1;
      if ( need_Pxx )
        Pxx = Pxx(1:psd_len) + [0; Pxx(Nfft:-1:psd_len+1); 0];
        if ( conf>0 )
          Vxx = Vxx(1:psd_len) + [0; Vxx(Nfft:-1:psd_len+1); 0];
        end
      end
      if ( need_Pxy )
        Pxy = Pxy(1:psd_len) + conj([0; Pxy(Nfft:-1:psd_len+1); 0]);
      end
      if ( need_Pyy )
        Pyy = Pyy(1:psd_len) + [0; Pyy(Nfft:-1:psd_len+1); 0];
      end
    else                    % one-sided, Nfft is odd
      psd_len = (Nfft+1)/2;
      if ( need_Pxx )
        Pxx = Pxx(1:psd_len) + [0; Pxx(Nfft:-1:psd_len+1)];
        if ( conf>0 )
          Vxx = Vxx(1:psd_len) + [0; Vxx(Nfft:-1:psd_len+1)];
        end
      end
      if ( need_Pxy )
        Pxy = Pxy(1:psd_len) + conj([0; Pxy(Nfft:-1:psd_len+1)]);
      end
      if ( need_Pyy )
        Pyy = Pyy(1:psd_len) + [0; Pyy(Nfft:-1:psd_len+1)];
      end
    end
  else                      % two-sided (and shifted)
    psd_len = Nfft;
  end
  % end MAIN CALCULATIONS
  %
  % SCALING AND OUTPUT
  % Put all results in matrix, one row per spectrum
  %   Pxx, Pxy, Pyy are sums of periodograms, so "n_ffts"
  %   in the scale factor converts them into averages
  spectra    = zeros(psd_len,n_results);
  spect_type = zeros(n_results,1);
  scale = n_ffts * seg_len * Fs * win_meansq;
  if ( do_power )
    spectra(:,do_power) = Pxx / scale;
    spect_type(do_power) = 1;
    if ( conf>0 )
      Vxx = [Pxx-Vxx Pxx+Vxx]/scale;
    end
  end
  if ( do_cross )
    spectra(:,do_cross) = Pxy / scale;
    spect_type(do_cross) = 2;
  end
  if ( do_trans )
    spectra(:,do_trans) = Pxy ./ Pxx;
    spect_type(do_trans) = 3;
  end
  if ( do_coher )
    % force coherence to be real
    spectra(:,do_coher) = real(Pxy .* conj(Pxy)) ./ Pxx ./ Pyy;
    spect_type(do_coher) = 4;
  end
  if ( do_ypower )
    spectra(:,do_ypower) = Pyy / scale;
    spect_type(do_ypower) = 5;
  end
  freq = [0:psd_len-1].' * ( Fs / Nfft );
  %
  % range='shift': Shift zero-frequency to the middle
  if ( range == 2 )
    len2 = fix((Nfft+1)/2);
    spectra = [ spectra(len2+1:Nfft,:); spectra(1:len2,:)];
    freq    = [ freq(len2+1:Nfft)-Fs; freq(1:len2)];
    if ( conf>0 )
      Vxx = [ Vxx(len2+1:Nfft,:); Vxx(1:len2,:)];
    end
  end
  %
  %  RETURN RESULTS or PLOT
  if ( nargout>=2 && conf>0 )
    varargout{2} = Vxx;
  end
  if ( nargout>=(2+(conf>0)) )
    % frequency is 2nd or 3rd return value,
    % depends on if 2nd is confidence interval
    varargout{2+(conf>0)} = freq;
  end
  if ( nargout>=1 )
    varargout{1} = spectra;
  else
    %
    % Plot the spectra if there are no return variables.
    plot_title=['power spectrum x ';
                'cross spectrum   ';
                'transfer function';
                'coherence        ';
                'power spectrum y ' ];
    for ii = 1: n_results
      if ( conf>0 && spect_type(ii)==1 )
        Vxxxx = Vxx;
      else
        Vxxxx = [];
      end
      if ( n_results > 1 )
        figure();
        end
      if ( plot_type == 1 )
        plot(freq,[abs(spectra(:,ii)) Vxxxx]);
      elseif ( plot_type == 2 )
        semilogx(freq,[abs(spectra(:,ii)) Vxxxx]);
      elseif ( plot_type == 3 )
        semilogy(freq,[abs(spectra(:,ii)) Vxxxx]);
      elseif ( plot_type == 4 )
        loglog(freq,[abs(spectra(:,ii)) Vxxxx]);
      elseif ( plot_type == 5 )  % db
        ylabel( 'amplitude (dB)' );
        plot(freq,[10*log10(abs(spectra(:,ii))) 10*log10(abs(Vxxxx))]);
      end
      title( char(plot_title(spect_type(ii),:)) );
      ylabel( 'amplitude' );
      % Plot phase of cross spectrum and transfer function
      if ( spect_type(ii)==2 || spect_type(ii)==3 )
        figure();
        if ( plot_type==2 || plot_type==4 )
          semilogx(freq,180/pi*angle(spectra(:,ii)));
        else
          plot(freq,180/pi*angle(spectra(:,ii)));
        end
        title( char(plot_title(spect_type(ii),:)) );
        ylabel( 'phase' );
      end
    end %for
  end 
end
end

%!demo
%! fflush(stdout);
%! rand('seed',2038014164);
%! a = [ 1.0 -1.6216505 1.1102795 -0.4621741 0.2075552 -0.018756746 ];
%! white = rand(1,16384);
%! signal = detrend(filter(0.70181,a,white));
%! % frequency shift by modulating with exp(j.omega.t) 
%! skewed = signal.*exp(2*pi*i*2/25*[1:16384]);
%! Fs = 25; % sampling frequency
%! hold off
%! pwelch([]);
%! pwelch(signal);
%! disp('Default settings: Fs=1Hz, overlap=0.5, no padding' )
%! input('Onesided power spectral density (real data). Press ENTER', 's' );
%! hold on
%! pwelch(skewed);
%! disp('Frequency-shifted complex data.  Twosided wrap-around spectrum.' );
%! input('Area is same as one-sided spectrum. Press ENTER', 's' );
%! pwelch(signal,'shift','semilogy');
%! input('Twosided, centred zero-frequency, lin-log plot. Press ENTER', 's' );
%! hold off
%! figure();
%! pwelch(skewed,[],[],[],Fs,'shift','semilogy');
%! input('Actual Fs=25 Hz. Note change of scales. Press ENTER', 's' );
%! pwelch(skewed,[],[],[],Fs,0.95,'shift','semilogy');
%! input('Spectral density with 95% confidence interval. Press ENTER', 's' );
%! pwelch('R12+');
%! pwelch(signal,'squared');
%! input('Spectral density with Matlab R12 defaults. Press ENTER', 's' );
%! figure();
%! pwelch([]);
%! pwelch(signal,3640,[],4096,2*pi,[],'no-strip');
%! input('Same spectrum with 95% confidence interval. Press ENTER', 's' );
%! figure();
%! pwelch(signal,[],[],[],2*pi,0.95,'no-strip');
%! input('95% confidence interval with native defaults. Press ENTER', 's' );
%! pwelch(signal,64,[],[],2*pi,'no-strip');
%! input('Only 32 frequency values in this spectrum. Press ENTER', 's' );
%! hold on
%! pwelch(signal,64,[],256,2*pi,'no-strip');
%! input('4:1 zero padding gives artificial smoothing. Press ENTER', 's' );
%! figure();
%! pwelch('psd');
%! pwelch(signal,'squared');
%! input('Just like Matlab spectrum.welch(...) defaults. Press ENTER', 's' );
%! hold off
%! pwelch({});
%! pwelch(white,signal,'trans','coher','short')
%! input('Transfer and coherence functions. Press ENTER', 's' );
%! disp('Use "close all" to remove plotting windows.' );

