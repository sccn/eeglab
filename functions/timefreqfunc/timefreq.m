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
%      'wavelet' = 0  -> Use FFTs (with constant window length) { Default } 
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
%      'detrend' = ['on'|'off'], Linearly detrend each data epoch  {'off'}
%
%    Optional FFT/DFT:
%      'tlimits'  = [min max] time limits in ms.
%      'winsize'  = If cycles==0 (FFT, see 'wavelet' input): data subwindow 
%                   length (fastest, 2^n<frames);
%                   if cycles >0: *longest* window length to use. This
%                   determines the lowest output frequency  {~frames/8}
%      'ntimesout' = Number of output times (int<frames-winsize). Enter a 
%                   negative value [-S] to subsample original time by S.
%      'timesout' = Enter an array to obtain spectral decomposition at 
%                   specific time values (note: algorithm find closest time 
%                   point in data and this might result in an unevenly spaced
%                   time array). Overwrite 'ntimesout'. {def: automatic}
%      'freqs'    = [min max] frequency limits. Default [minfreq srate/2], 
%                   minfreq being determined by the number of data points, 
%                   cycles and sampling frequency. Enter a single value
%                   to compute spectral decompisition at a single frequency
%                   (note: for FFT the closest frequency will be estimated).
%                   For wavelet, reducing the max frequency reduce
%                   the computation load.
%      'padratio' = FFTlength/winsize (2^k)                     {def: 2}
%                   Multiplies the number of output frequencies by
%                   dividing their spacing. When cycles==0, frequency
%                   spacing is (low_frequency/padratio).
%      'nfreqs'   = number of output frequencies. For FFT, closest computed
%                   frequency will be returned. Overwrite 'padratio' effects
%                   for wavelets. Default: use 'padratio'.
%      'freqscale' = ['log'|'linear'] frequency scale. Default is 'linear'.
%                    Note that for obtaining 'log' spaced freqs using FFT, 
%                    closest correspondant frequencies in the 'linear' space 
%                    are returned.
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
% Note: it is not advised to use a FFT decomposition in a log scale. Output
%       value are accurate but plotting might not be because of the non-uniform
%       frequency output values in log-space. If you have to do it, use a 
%       padratio as large as possible, or interpolate time-freq image at
%       exact log scale values before plotting.
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
% Revision 1.43  2004/02/27 19:07:07  arno
% fixing fft limits
%
% Revision 1.42  2004/02/24 18:46:33  arno
% debug FFT log ...
%
% Revision 1.41  2004/02/17 19:56:05  arno
% debug frequency assigment for FFT
%
% Revision 1.40  2003/12/10 18:15:11  arno
% msg
%
% Revision 1.39  2003/11/24 20:58:50  arno
% debug last modif
%
% Revision 1.38  2003/11/22 02:28:15  arno
% msg for removing time points
%
% Revision 1.37  2003/11/12 02:49:00  arno
% adding warnings, reprograming freqs
%
% Revision 1.36  2003/10/30 18:46:06  arno
% time limits in ms
%
% Revision 1.35  2003/10/29 01:43:27  arno
% typo in error msg
%
% Revision 1.34  2003/10/29 01:22:06  arno
% remve debug msg
%
% Revision 1.33  2003/10/29 01:20:21  arno
% more error checking
%
% Revision 1.32  2003/10/29 00:13:57  arno
% processing a single frequency
%
% Revision 1.31  2003/10/22 18:07:37  arno
% if only 1 input for g.freqs
%
% Revision 1.30  2003/10/22 17:48:17  arno
% debuging default spacing between frequencies
% ,
%
% Revision 1.29  2003/10/22 17:06:54  arno
% winsize length
%
% Revision 1.28  2003/07/09 21:29:57  arno
% implementing ntimesout
%
% Revision 1.27  2003/06/27 15:29:45  arno
% removing deprecated arguments
%
% Revision 1.26  2003/06/27 01:00:26  arno
% updating timeout help
%
% Revision 1.25  2003/06/19 00:24:01  arno
% adding error message if negative number of time points
%
% Revision 1.24  2003/05/29 15:14:30  arno
% allowing vector for timesout
%
% Revision 1.23  2003/05/29 14:59:07  arno
% nothing
%
% Revision 1.22  2003/05/22 01:18:56  arno
% header typo
%
% Revision 1.21  2003/05/02 17:57:42  arno
% removing maxfreq
%
% Revision 1.20  2003/05/02 17:34:12  arno
% implementing ITC option for very big arrays
%
% Revision 1.19  2003/04/30 00:07:25  arno
% debug if window size changes
%
% Revision 1.18  2003/04/29 21:53:35  arno
% rounding trials
%
% Revision 1.17  2003/04/29 18:41:27  arno
% debug computation
% .,
%
% Revision 1.16  2003/04/29 01:08:44  arno
% new version which can generate wavelets at any frequencies
%
% Revision 1.15  2003/01/10 01:26:59  arno
% nothing
%
% Revision 1.14  2003/01/08 23:32:04  arno
% typo in linear coherence formula
%
% Revision 1.13  2003/01/06 17:52:38  arno
% nothing
%
% Revision 1.12  2003/01/05 18:55:12  arno
% testing
%
% Revision 1.11  2003/01/02 20:01:39  arno
% removing the complex transpose bug
%
% Revision 1.10  2002/10/18 14:43:35  arno
% header
%
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

function [tmpall, freqs, timesout, itcvals] = timefreq(data, srate, varargin)
%	nb_points, timesout, naccu, baselength, baseboot, boottype, alpha, rboot);
	
if nargin < 2
	help timefreq;
	return;
end;

[frame trials]= size(data);
g = finputcheck(varargin, ...
                { 'ntimesout'     'integer'  []                       200; ...
				  'timesout'      'real'     []                       []; ...
				  'winsize'       'integer'  [0 Inf]                  max(pow2(nextpow2(frame)-3),4); ...
                  'tlimits'       'real'     []                       [0 frame/srate*1000]; ...
                  'detrend'       'string'   {'on' 'off'}              'off'; ...
				  'freqs'         'real'     [0 Inf]                  [0 srate/2]; ...
				  'nfreqs'        'integer'  [0 Inf]                  []; ...
				  'freqscale'     'string'   { 'linear' 'log' }       'linear'; ...
				  'wavelet'       'real'     [0 Inf]                   0; ...
				  'padratio'      'integer'  [1 Inf]                   2; ...
				  'itctype'       'string'   {'phasecoher' 'phasecoher2' 'coher'}  'phasecoher'; ...
				  'subitc'        'string'   {'on' 'off'}              'off'	});
if isstr(g), error(g); end;

% checkin parameters
% ------------------
g.cycles = g.wavelet(1);
if length(g.wavelet) >= 2
	g.cyclesfact = g.wavelet(2);
else 
	g.cyclesfact = 1; % default is 1 (no decrease)
end;
if g.cyclesfact < 0 | g.cyclesfact > 1
    error('Wavelet cycle factor must be comprised between 0 and 1');
end;
if (g.cycles == 0 & pow2(nextpow2(g.winsize)) ~= g.winsize)
   error('Value of winsize must be an integer power of two [1,2,4,8,16,...]');
elseif (g.winsize > frame)
   error('Value of winsize must be less than frame length.');
end     
if (pow2(nextpow2(g.padratio)) ~= g.padratio)
   error('Value of padratio must be an integer power of two [1,2,4,8,16,...]');
end

% finding frequency limits
% ------------------------
if g.cycles ~= 0 & g.freqs(1) == 0, g.freqs(1) = srate*g.cycles/g.winsize; end;

% finding frequencies
% -------------------
if length(g.freqs) == 2
    
    % min and max
    % -----------
    if g.freqs(1) == 0 & g.cycles ~= 0
        g.freqs(1) = srate*g.cycles/g.winsize;
    end;

    % default number of freqs using padratio
    % --------------------------------------
    if isempty(g.nfreqs)
        g.nfreqs = g.winsize/2*g.padratio+1;
        % adjust nfreqs depending on frequency range
        tmpfreqs = linspace(0, srate/2, g.nfreqs); 
        tmpfreqs = tmpfreqs(2:end);  % remove DC (match the output of PSD)
        
        % adjust limits for FFT (only linear scale)
        if g.cycles == 0 & ~strcmpi(g.freqscale, 'log')
            if ~any(tmpfreqs == g.freqs(1))
                [tmp minind] = min(abs(tmpfreqs-g.freqs(1)));
                g.freqs(1)   = tmpfreqs(minind);
                fprintf('Adjust min freq. to %3.2f Hz to match FFT output frequencies\n', g.freqs(1));
            end;
            if ~any(tmpfreqs == g.freqs(2))
                [tmp minind] = min(abs(tmpfreqs-g.freqs(2)));
                g.freqs(2)   = tmpfreqs(minind);
                fprintf('Adjust max freq. to %3.2f Hz to match FFT output frequencies\n', g.freqs(2));
            end;
        end;
        
        % find number of frequencies
        % --------------------------
        g.nfreqs = length(tmpfreqs( intersect( find(tmpfreqs >= g.freqs(1)), find(tmpfreqs <= g.freqs(2)))));
        if g.freqs(1)==g.freqs(2), g.nfreqs = 1; end;
    end;

    % find closest freqs for FFT
    % --------------------------
    if strcmpi(g.freqscale, 'log')
        g.freqs = linspace(log(g.freqs(1)), log(g.freqs(end)), g.nfreqs);
        g.freqs = exp(g.freqs);
    else
        g.freqs = linspace(g.freqs(1), g.freqs(2), g.nfreqs); % this should be OK for FFT
                                                              % because of the limit adjustment
    end;
end;
g.nfreqs = length(g.freqs);


% function for time freq initialisation
% -------------------------------------
if (g.cycles == 0) %%%%%%%%%%%%%% constant window-length FFTs %%%%%%%%%%%%%%%%
    freqs = linspace(0, srate/2, g.winsize*g.padratio/2+1);
    freqs = freqs(2:end); % remove DC (match the output of PSD)
                          %srate/g.winsize*[1:2/g.padratio:g.winsize]/2
    g.win   = hanning(g.winsize);
else % %%%%%%%%%%%%%%%%%% Constant-Q (wavelet) DFTs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %freqs = srate*g.cycles/g.winsize*[2:2/g.padratio:g.winsize]/2;
     %g.win = dftfilt(g.winsize,g.freqs(2)/srate,g.cycles,g.padratio,g.cyclesfact);

    freqs = g.freqs;
    if g.cyclesfact ~= 1, 
        g.cycles = [ g.cycles g.cycles*g.freqs(end)/g.freqs(1)*(1-g.cyclesfact)]; 
    end;
    g.win    = dftfilt2(g.freqs,g.cycles,srate, g.freqscale);
    
    g.winsize = 0;
    for index = 1:length(g.win)
        g.winsize = max(g.winsize,length(g.win{index}));
    end;
end;

% compute time vector
% -------------------
[ g.timesout g.indexout ] = gettimes(frame, g.tlimits, g.timesout, g.winsize, g.ntimesout);
tmpall      = repmat(nan,[length(freqs) length(g.timesout) trials]);

% -------------------------------
% compute time freq decomposition
% -------------------------------
fprintf('The window size used is %d samples (%g ms) wide.\n',g.winsize, 1000/srate*g.winsize);
fprintf('Estimating %d %s-spaced frequencies from %2.1f Hz to %3.1f Hz.\n', length(g.freqs), ...
        fastif(strcmpi(g.freqscale, 'log'), 'log', 'linear'), g.freqs(1), g.freqs(end));

if g.cycles(1) == 0
    fprintf('Processing trial (of %d):',trials);
    for trial = 1:trials
        if rem(trial,10) == 0,  fprintf(' %d',trial); end
        if rem(trial,120) == 0, fprintf('\n'); end
        for index = 1:length(g.indexout)
            tmpX = data([-g.winsize/2+1:g.winsize/2]+g.indexout(index)+(trial-1)*frame); % 1 point imprecision
            
            tmpX = tmpX - mean(tmpX);
            if strcmpi(g.detrend, 'on'),
                tmpX = detrend(tmpX); 
            end;
            
            tmpX = g.win .* tmpX(:);
            tmpX = fft(tmpX,g.padratio*g.winsize);
            tmpX = tmpX(2:g.padratio*g.winsize/2+1);
            tmpall(:,index, trial) = tmpX(:);            
        end;
    end;
else % wavelet
    % prepare wavelet filters
    % -----------------------
    for index = 1:length(g.win)
        g.win{index} = transpose(repmat(g.win{index}, [trials 1]));
    end;
    
    % apply filters
    % -------------
    fprintf('Processing time point (of %d):',length(g.timesout));
    for index = 1:length(g.indexout)
        if rem(index,10) == 0,  fprintf(' %d',index); end
        if rem(index,120) == 0, fprintf('\n'); end
        for freqind = 1:length(g.win)
            wav = g.win{freqind}; sizewav = size(wav,1)-1;
            %g.indexout(index), size(wav,1), g.freqs(freqind)
            tmpX = data([-sizewav/2:sizewav/2]+g.indexout(index),:);
            
            tmpX = tmpX - ones(size(tmpX,1),1)*mean(tmpX);
            if strcmpi(g.detrend, 'on'),
                for trial = 1:trials
                    tmpX(:,trial) = detrend(tmpX(:,trial)); 
                end;
            end;
            
            tmpX = sum(wav .* tmpX);
            tmpall( freqind, index, :) = tmpX;            
        end;
    end;
end;    
fprintf('\n');

% compute and subtract ITC
% ------------------------
if nargout > 3 | strcmpi(g.subitc, 'on')
	itcvals = tfitc(tmpall, g.itctype);
end;
if strcmpi(g.subitc, 'on')
    %a = gcf; figure; imagesc(abs(itcvals)); cbar; figure(a);
	tmpall = (tmpall - abs(tmpall) .* repmat(itcvals, [1 1 trials])) ./ abs(tmpall);
end;

% find closest output frequencies
% -------------------------------
if length(g.freqs) ~= length(freqs) | any(g.freqs ~= freqs)
    allindices = zeros(1,length(g.freqs));
    for index = 1:length(g.freqs)
        [dum ind] = min(abs(freqs-g.freqs(index)));
        allindices(index) = ind;
    end;
    fprintf('finding closest frequencies: %d freqs removed\n', length(freqs)-length(allindices));
    freqs = freqs(allindices);
    tmpall = tmpall(allindices,:,:);
    if nargout > 3 | strcmpi(g.subitc, 'on')
        itcvals = itcvals(allindices,:,:);
    end;
end;    

timesout = g.timesout;

%figure; imagesc(abs(sum(itcvals,3))); cbar;
return;


% function for itc
% ----------------
function [itcvals] = tfitc(tfdecomp, itctype);
% first dimension are trials
switch itctype
 case 'coher',
  try,
      itcvals = sum(tfdecomp,3) ./ sqrt(sum(tfdecomp .* conj(tfdecomp),3) * size(tfdecomp,3));
  catch, % scan rows if out of memory
      for index =1:size(tfdecomp,1)
          itcvals(index,:) = sum(tfdecomp(index,:,:),3) ./ sqrt(sum(tfdecomp(index,:,:) .* conj(tfdecomp(index,:,:)),3) * size(tfdecomp,3));
      end;
  end;
 case 'phasecoher2',
  try,
      itcvals = sum(tfdecomp,3) ./ sum(sqrt(tfdecomp .* conj(tfdecomp)),3);
  catch, % scan rows if out of memory
      for index =1:size(tfdecomp,1)
          itcvals(index,:) = sum(tfdecomp(index,:,:),3) ./ sum(sqrt(tfdecomp(index,:,:) .* conj(tfdecomp(index,:,:))),3);
      end;
  end;
 case 'phasecoher',
  try,
      itcvals = sum(tfdecomp ./ sqrt(tfdecomp .* conj(tfdecomp)) ,3) / size(tfdecomp,3);
  catch, % scan rows if out of memory
      for index =1:size(tfdecomp,1)
          itcvals(index,:) = sum(tfdecomp(index,:,:) ./ sqrt(tfdecomp(index,:,:) .* conj(tfdecomp(index,:,:))) ,3) / size(tfdecomp,3);
      end;
  end;
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

% get time points
% ---------------
function [ timevals, timeindices ] = gettimes(frames, tlimits, timevar, winsize, ntimevar);
    timevect = linspace(tlimits(1), tlimits(2), frames);
    srate    = 1000*(frames-1)/(tlimits(2)-tlimits(1));
    
    if isempty(timevar) % no pre-defined time points 
        if ntimevar(1) > 0
            % generate linearly space vector
            % ------------------------------
            if (ntimevar > frames-winsize)
                ntimevar = frames-winsize;
                if ntimevar < 0
                    error('Not enough data points, reduce the window size or lowest frequency');
                end;
                disp(['Value of ''timesout'' must be <= frame-winsize, ''timesout'' adjusted to ' int2str(ntimevar) ]);
            end
            npoints = ntimevar(1);
            wintime = 500*winsize/srate;
            timevals = linspace(tlimits(1)+wintime, tlimits(2)-wintime, npoints);
            fprintf('Generating %d time points (%1.1f to %1.1f ms)\n', npoints, min(timevals), max(timevals));
        else            
            % subsample data
            % --------------
            nsub     = -ntimevar(1);
            timeindices = [ceil(winsize/2+nsub/2):nsub:length(timevect)-ceil(winsize/2)-1];
            timevals    = timevect( timeindices );
            fprintf('Subsampling by %d (%1.1f to %1.1f ms)\n', nsub, min(timevals), max(timevals));
        end;
    else
        timevals = timevar;
        % check boundaries
        % ----------------
        wintime = 500*winsize/srate;
        tmpind  = find( (timevals >= tlimits(1)+wintime) & (timevals <= tlimits(2)-wintime) );
        if  length(timevals) ~= length(tmpind)
            fprintf('Warning: %d out of %d time values were removed (now %3.2f to %3.2f ms) so the lowest\n', ...
                    length(timevals)-length(tmpind), length(timevals), timevals(tmpind(1)), timevals(tmpind(end)));
            fprintf('         frequency could be computed with the requested accuracy\n');
        end;        
        timevals = timevals(tmpind);
    end;
    
    % find closet points in data
    % --------------------------
    oldtimevals = timevals;
    for index = 1:length(timevals)
        [dum ind] = min(abs(timevect-timevals(index)));
        timeindices(index) = ind;
        timevals(index)    = timevect(ind);
    end;
    if length(timevals) < length(unique(timevals))
        disp('Warning: duplicate times, reduce the number of output times');
    end;
    if all(oldtimevals == timevals)
        disp('Debug msg: Time value unchanged by finding closest in data');        
    else 
        %disp('Debug msg: Time value updated by finding closest points in data');
    end;    

% DEPRECATED, FOR C INTERFACE
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

