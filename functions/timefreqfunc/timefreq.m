% timefreq() - compute time/frequency decomposition of data trials.
%
% Usage:
%     >> [accarray, rsignif] = timefreq(arg1, arg2, formula, varargin ...);
%
% Inputs:
%    arg1    - [array] 2-D or 3-D array of value
%    arg2    - [array] 2-D or 3-D array of value
%    formula - [string] formula to compute a given measure. Takes arguments
%              'arg1', 'arg2' as inputs and 'res' (by default) as output. i.e.
%              'res = res + arg1 .* conj(arg2)'
%
% Optional inputs:
%   'boottype '   - ['first'|'second'|'both'] 'fist' accumulate in the first 
%                 dimension only (shuffling the first dimension). 'second'
%                 shuffle the second dimension but keep the fist constant. 
%                 'both' shuffle the two dimension. Default is 'both'.
%   'alpha'       - [real] significance (between 0 and 1). Default is 0.05.
%   'naccu'       - [integer] number of exemplar to accumulate. Default is 200.
%   'basevect'    - [integer vector] time vector indices for baseline. Default 
%                 is all time points.
%   'accarray'    - accumulation array (from previous calls). Allow to compute
%                 only the 'rsignif' output.
%   'formulainit' - [string] for initializing variable. i.e. 'res = zeros(10,40);'
%                 Default is initialization of the res variable to
%                 (size(arg1,3) x naccu)
%   'formulapost' - [string] after the accumulation. i.e. 'res = res /10;'
%                 default is none.
%   'formulaout'  - [string] name of the computed variable. Default is 'res'.
%
% Outputs: 
%    diffres  - differential array computed on the non-shuffled data
%    accdres  - result for shuffled data
%    res1     - result for first condition
%    res2     - result for second condition
%
% Authors: Arnaud Delorme, Lars & Scott Makeig
%          CNL/Salk Institute 1998-2001; SCCN/INC/UCSD, La Jolla, 2002-
%
% See also: timef()

% Copyright (C) 8/1/98  Arnaud Delorme, Sigurd Enghoff & Scott Makeig, SCCN/INC/UCSD
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: not supported by cvs2svn $
% Revision 1.2  2002/10/02 00:36:08  arno
% update condstat, debug
%
% Revision 1.1  2002/10/01 16:09:52  arno
% Initial revision
%
% Revision 1.1  2002/09/24 23:28:02  arno
% Initial revision
%
% function for bootstrap initialisation
% -------------------------------------

function [tmpall, freqs, times, itcvals] = timefreq(data, srate, varargin)
%	nb_points, timesout, naccu, baselength, baseboot, boottype, alpha, rboot);
	
if nargin < 2
	help timefreq;
	return;
end;

[frame trials]= size(data);
g = finputcheck(varargin, ...
                { 'timesout'      'integer'  [0 Inf]                  200; ...
				  'winsize'       'integer'  [0 Inf]                  max(pow2(nextpow2(frame)-3),4); ...
                  'tlimits'       'real'     []                       [0 frame/srate]; ...
                  'basevect'      'integer'  [0 Inf]                  []; ...
                  'detrend'       'string'   {'on' 'off'}              'off'; ...
				  'maxfreq'       'real'     [0 Inf]                  (srate/2); ...
				  'freq'          'real'     [0 Inf]                  []; ...
				  'trial'         'integer'  [0 Inf]                  []; ...
				  'wavelet'       'real'     [0 Inf]                   0; ...
				  'padratio'      'integer'  [1 Inf]                   2; ...
				  'mtaper'        'real'     []                        []; ...
				  'itctype'       'string'   {'phasecoher' 'phasecoher2' 'coher'}  'phasecoher'; ...
				  'subitc'        'string'   {'on' 'off'}              'off'	});

g.cycles = g.wavelet(1);
if length(g.wavelet) >= 2
	g.cyclesfact = g.wavelet(2);
else 
	g.cyclesfact = 0.5;
end;
if (g.cycles == 0 & pow2(nextpow2(g.winsize)) ~= g.winsize)
   error('Value of winsize must be an integer power of two [1,2,4,8,16,...]');
elseif (g.winsize > frame)
   error('Value of winsize must be less than frame length.');
end 
if (g.timesout > frame-g.winsize)
   g.timesout = g.frame-g.winsize;
   disp(['Value of timesout must be <= frame-winsize, timeout adjusted to ' int2str(g.timesout) ]);
end
if (pow2(nextpow2(g.padratio)) ~= g.padratio)
   error('Value of padratio must be an integer power of two [1,2,4,8,16,...]');
end
if (g.maxfreq > srate/2)
   fprintf('Warning: value of g.maxfreq greater that Nyquist rate\n\n');
end

% compute time limits
% -------------------
wintime = 500*g.winsize/srate;
times = [g.tlimits(1)+wintime:(g.tlimits(2)-g.tlimits(1)-2*wintime)/(g.timesout-1):g.tlimits(2)-wintime];

% function for time freq initialisation
% -------------------------------------
g.stp       = (frame-g.winsize)/(g.timesout-1);
if (g.cycles == 0) %%%%%%%%%%%%%% constant window-length FFTs %%%%%%%%%%%%%%%%
   freqs = srate/g.winsize*[1:2/g.padratio:g.winsize]/2;
   g.win   = hanning(g.winsize);
   g.nb_points = g.padratio*g.winsize/2;   
else % %%%%%%%%%%%%%%%%%% Constant-Q (wavelet) DFTs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   freqs = srate*g.cycles/g.winsize*[2:2/g.padratio:g.winsize]/2;
   g.win = dftfilt(g.winsize,g.maxfreq/srate,g.cycles,g.padratio,g.cyclesfact);
   g.nb_points = size(g.win,2);
end;
%tmpall      = repmat(nan,[trials g.timesout g.nb_points]);
tmpall      = repmat(nan,[g.nb_points g.timesout trials]);

% ------------------------------------
% compute time freq decomposition
% ------------------------------------
for trial = 1:trials
	if rem(trial,10) == 0,  fprintf(' %d',trial); end
	if rem(trial,120) == 0, fprintf('\n'); end
	for index = 1:g.timesout
		tmpX = data([1:g.winsize]+floor((index-1)*g.stp)+(trial-1)*frame);
		
		tmpX = tmpX - mean(tmpX);
		if strcmpi(g.detrend, 'on'),
			tmpX = detrend(tmpX); 
		end;
            
		if g.cycles == 0 % use FFTs
			tmpX = g.win .* tmpX(:);
			tmpX = fft(tmpX,g.padratio*g.winsize);
			tmpX = tmpX(2:g.padratio*g.winsize/2+1);
		else 
			tmpX = g.win' * tmpX(:);
		end
		tmpall(:,index, trial) = tmpX(:);

		%if g.ITCdone
        %    tmpX = (tmpX - abs(tmpX) .* g.ITC(:,index)) ./ abs(tmpX);
		%end;
   end;
end;

% compute and subtract ITC
% ------------------------
if nargout > 3 | strcmpi(g.subitc, 'on')
	itcvals = tfitc(tmpall, g.itctype);
end;
if strcmpi(g.subitc, 'on')
	%itcvals = transpose(itcvals); % do not use ' otherwise conjugate
	
	itcvalsub = repmat(shiftdim(itcvals, -1), [trials 1 1]);
	tmpall = (tmpall - abs(tmpall) .* itcvalsub) ./ abs(tmpall);
end;
return;


% function for itc
% ----------------
function [itcvals] = tfitc(tfdecomp, itctype);
% first dimension are trials
switch itctype
   case 'coher',
      itcvals = sum(tfdecomp,3) ./ sqrt(sum(tfdecomp .* conj(tfdecomp),3)) / size(tfdecomp,3);
	  %g.ITC(:,times)      = g.ITC(:,times) + g.tmpalltimes; % complex coher.
      %g.ITCcumul(:,times) = g.ITCcumul(:,times)+abs(g.tmpalltimes).^2;
      %case 'coher',       g.ITC = g.ITC ./ sqrt(trials * g.ITCcumul);
   case 'phasecoher2',
      itcvals = sum(tfdecomp,3) ./ sum(sqrt(tfdecomp .* conj(tfdecomp)),3);
      %g.ITC(:,times)      = g.ITC(:,times) + g.tmpalltimes; % complex coher.
      %g.ITCcumul(:,times) = g.ITCcumul(:,times)+abs(g.tmpalltimes);
      %case 'phasecoher2', g.ITC = g.ITC ./ g.ITCcumul;
   case 'phasecoher',
      itcvals = sum(tfdecomp ./ sqrt(tfdecomp .* conj(tfdecomp)) ,3) / size(tfdecomp,3);
      %g.ITC(:,times)      = g.ITC(:,times) + g.tmpalltimes ./ abs(g.tmpalltimes); % complex coher.
      %case 'phasecoher',  g.ITC = g.ITC / trials; % complex coher.
end % ~any(isnan())
return;

function w = hanning(n)
if ~rem(n,2)
   w = .5*(1 - cos(2*pi*(1:n/2)'/(n+1)));
   w = [w; w(end:-1:1)];
else
   w = .5*(1 - cos(2*pi*(1:(n+1)/2)'/(n+1)));
   w = [w; w(end-1:-1:1)];
end

function nofunction()
   % C PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   filename = [ 'tmpcrossf' num2str(round(rand(1)*1000)) ];
   f = fopen([ filename '.in'], 'w');
   fwrite(f, tmpsaveall, 'int32');
   fwrite(f, g.detret, 'int32');
   fwrite(f, g.srate, 'int32');
   fwrite(f, g.maxfreq, 'int32');
   fwrite(f, g.padratio, 'int32');
   fwrite(f, g.cycles, 'int32');
   fwrite(f, g.winsize, 'int32');
   fwrite(f, g.timesout, 'int32');
   fwrite(f, g.subitc, 'int32');
   fwrite(f, g.type, 'int32');
   fwrite(f, trials, 'int32');
   fwrite(f, g.naccu, 'int32');
   fwrite(f, length(X), 'int32');
   fwrite(f, X, 'double');
   fwrite(f, Y, 'double');
   fclose(f);
   
   command = [ '!cppcrosff ' filename '.in ' filename '.out' ];
   eval(command);
   
   f = fopen([ filename '.out'], 'r');
   size1 = fread(f, 'int32', 1);
   size2 = fread(f, 'int32', 1);
   Rreal = fread(f, 'double', [size1 size2]);
   Rimg  = fread(f, 'double', [size1 size2]);
   Coher.R = Rreal + j*Rimg;
   Boot.Coherboot.R = [];
   Boot.Rsignif = [];

