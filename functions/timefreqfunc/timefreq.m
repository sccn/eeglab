% timefreq() - compute time/frequency decomposition of data trials.
%
% Usage:
%     >> [tf, freqs, times] = timefreq(data, srate);
%     >> [tf, freqs, times, itcvals] = timefreq(data, srate, 'key1', 'val1', ...
%                                      'key2', 'val2' ...)
%
% Inputs:
%       data    = [float array] 2-D data array of size (times,trials)
%       srate   = sampling rate
%
% Optional inputs:
%       wavelet = 0  -> Use FFTs (with constant window length) { Default } 
%               = >0 -> Number of cycles in each analysis wavelet 
%               = [cycles expfactor] -> if 0 < expfactor < 1,  the number 
%                 of wavelet cycles expands with frequency from cycles
%                 If expfactor = 1, no expansion; if = 0, constant
%                 window length (as in FFT)            {default wavelet: 0}
%
%    Optional Coherence Type for ITC:
%      'type'   = ['coher'|'phasecoher'] Compute either linear coherence
%                 ('coher') or phase coherence ('phasecoher') also known
%                 as phase coupling factor' {default: 'phasecoher'}.
%      'subitc' = ['on'|'off'] subtract stimulus locked Inter-Trial Coherence 
%                 (ITC) from x and y. This computes the  'intrinsic' coherence
%                 x and y not arising from common synchronization to 
%                 experimental events. See notes. {default: 'off'}
%
%    Optional Detrend:
%      'detrep' = ['on'|'off'], Linearly detrend data across trials  {'off'}
%
%    Optional FFT/DFT:
%      'winsize'  = If cycles==0 (FFT, see 'wavelet' input): data subwindow 
%                   length (fastest, 2^n<frames);
%                   if cycles >0: *longest* window length to use. This
%                   determines the lowest output frequency  {~frames/8}
%      'timesout' = Number of output times (int<frames-winsize) {def: 200}
%      'padratio' = FFTlength/winsize (2^k)                     {def: 2}
%                    Multiplies the number of output frequencies by
%                    dividing their spacing. When cycles==0, frequency
%                    spacing is (low_frequency/padratio).
%      'maxfreq'  = Maximum frequency (Hz) to plot (& output if cycles>0) 
%                    If wavelet==0, all FFT frequencies are output.
%                    For wavelet, reducing the max frequency reduce
%                    the computation load {def: 50}
%
% Outputs: 
%       tf      - time frequency array for all trials (freqs, times, trials)
%       freqs   - vector of computed frequencies (Hz)
%       times   - vector of computed time points (ms)
%       itcvals - time frequency "average" for all trials (freqs, times).
%                 In the coherence case, it is the real mean of the time
%                 frequency decomposition, but in the phase coherence case
%                 (see 'type' input'), this is the mean of the normalized
%                 spectral estimate.
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
% Revision 1.9  2002/10/16 01:21:16  arno
% removing debug plots
%
% Revision 1.8  2002/10/16 00:18:34  arno
% debuging with cooper
%
% Revision 1.7  2002/10/15 23:47:29  arno
% subitc options
%
% Revision 1.6  2002/10/15 22:58:59  arno
% nothing
%
% Revision 1.5  2002/10/15 18:29:46  arno
% changing default wavelet factor
%
% Revision 1.4  2002/10/15 00:08:50  arno
% g.frame -> frame
%
% Revision 1.3  2002/10/09 00:03:20  arno
% debugging
%
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

% checkin parameters
% ------------------
g.cycles = g.wavelet(1);
if length(g.wavelet) >= 2
	g.cyclesfact = g.wavelet(2);
else 
	g.cyclesfact = 1; % default is 1 (no decrease)
end;
if (g.cycles == 0 & pow2(nextpow2(g.winsize)) ~= g.winsize)
   error('Value of winsize must be an integer power of two [1,2,4,8,16,...]');
elseif (g.winsize > frame)
   error('Value of winsize must be less than frame length.');
end 
if (g.timesout > frame-g.winsize)
   g.timesout = frame-g.winsize;
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
    %a = gcf; figure; imagesc(abs(itcvals)); cbar; figure(a);
	tmpall = (tmpall - abs(tmpall) .* repmat(itcvals, [1 1 trials])) ./ abs(tmpall);
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

