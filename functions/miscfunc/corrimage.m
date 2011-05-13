% corrimage() - compute correlation image between an event and amplitude
%               and phase in the time-frequency domain.
%
% Usage: 
%     corrimage(data, sortvar, times);
%     [times, freqs, alpha, sigout, limits ] = corrimage(data, sortvar, ...
%                                    times, 'key1', val1, 'key2', val2, ...)
%
% Inputs:
%   data    - [float array] data of size (times,trials)
%   sortvar - [float array] sorting variable of size (trials)
%   times   - [float array] time indices for every data point. Size of
%             (times). Note: same input as the times vector for erpimage.
%
% Optional input:
%   'mode'     - ['phase'|'amp'] compute correlation of event values with 
%                phase ('phase') or amplitude ('amp') of signal at the 
%                given time-frequency points. Default is 'amp'.
%   'freqs'    - [min nfreqs max] build a frequency vector of size nfreqs 
%                with frequency spaced in a log scale. Then compute 
%                correlation at these frequency value.
%                Default is 50 points in a log scale between 2.5Hz to 50Hz.
%   'times'    - [float vector] compute correlation at these time value (ms).
%                (uses closest times points in 'times' input and return
%                them in the times output).
%                Default is 100 steps between min time+5% of time range and
%                max time-5%. Enter only 2 values [N X] to generate N 
%                time points and trim by X %. If N is negative, uses it as
%                a subsampling factor [-3 5] trims times by 5% and subsample
%                by 3 (by subsampling one obtains a regularly spaced times).
%   'trim'     - [low high] low and high percentile of sorted sortvar values
%                to retain. i.e. [5 95] remove the 5 lowest and highest 
%                percentile of sortvar values (and associated data) before
%                computing statistics. Default is [0 100].
%   'align'    - [float] same as 'align' parameter of erpimage(). This 
%                parameter is used to contrain the 'times' parameter so
%                correlation with data trials containing 0-values (as a
%                result of data alignment) are avoided: computing these
%                correlations would produce spurious significant results.
%                Default is no alignment.
%   'method'   - ['erpimage'|timefreq'] use either the erpimage() function
%                of the timefreq() function to compute spectral decomposition.
%                Default is 'timefreq' for speed reasons (note that both 
%                methods should return the same results).
%   'erpout'   - [min max] regress out ERP using the selected time-window [min
%                max] in ms (the ERP is subtracted from the whole time period 
%                but only regressed out in the selected time window).
%   'triallimit' - [integer array] specify trial boundaries for subjects. For 
%                instance [1 200 400] indicates 2 subjects, trials 1 to 199
%                for subject 1 and trials 200 to 399 for subject 2. This is 
%                currently only used for regressing erp out.
%
% Processing options:
%   'erpimopt' - [cell array] erpimage additional options (number of cycle ...).
%   'tfopt'    - [cell array] timefreq additional options (number of cycle ...).
% 
% Plotting options:
%   'plotvals' - [cell array of output values] enter output values here
%                { times freqs alpha sigout} to replot them.
%   'nofig'    - ['on'|'off'] do not create figure.
%   'cbar'     - ['on'|'off'] plot color bar. Default: 'on'.
%   'smooth'   - ['on'|'off'] smooth significance array. Default: 'on'.
%   'plot'     - ['no'|'alpha'|'sigout'|'sigoutm'|'sigoutm2'] plot -10*log(alpha)
%                values ('alpha'); output signal (slope or ITC) ('sigout'), 
%                output signal masked by significance ('sigoutm') or the last
%                2 option combined ('sigoutm2'). 'no' prevent the function
%                from plotting. In addition, see pmask. Default is 'sigoutmasked'.
%   'pmask'    - [real] maximum p value to show in plot. Default is 0.00001
%                (0.001 taking into account multiple comparisons (100)). Enter 
%                0.9XXX to get the higher tail or a negative value (e.e., -0.001
%                to get both tails).
%   'vert'     - [real array] time vector for vertivcal lines.
%   'limits'   - [min max] plotting limits. 
%
% Outputs:
%   times    - [float vector] vector of times (ms)
%   freqs    - [float vector] vector of frequencies (Hz)
%   alpha    - [float array] array (freqs,times) of p-values.
%   sigout   - [float array] array (freqs,times) of signal out (coherence
%              for phase and slope for amplitude).
%
% Important note: the 'timefreq' method may truncate the time window to
%                 compute the spectral decomposition at the lowest freq.
%
% Author: Arnaud Delorme & Scott Makeig, SCCN UCSD, 
%         and CNL Salk Institute, 18 April 2003
%
% See also: erpimage()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2002 Arnaud Delorme
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

function [time, freq, alpha, sigout, limits, tf, sortvar] = corrimage(data, sortvar, timevect, varargin)

if nargin < 3
    help corrimage;
    return;
end;
    
if isempty(timevect), timevect = [1 2]; end;

% check inputs
% ------------
g = finputcheck(varargin, { 'freqs'    'real'   [0 Inf]    [2.5 50 50];
                            'times'    'real'   []         [100 5]; % see function at the end
                            'mode'     'string' { 'phase','amp' }  'amp'; 
                            'vert'     'real'   []         [];
                            'align'    { 'real','cell' }   []         []; 
                            'plotvals' 'cell'   []         {}; 
                            'pmask'    'real'   []         0.00001; 
                            'triallimit' 'integer'   []       []; 
                            'trim'     'real'   [0 100]    [0 100]; 
                            'limits'   'real'   []         [];
                            'method'   'string' { 'erpimage','timefreq' }        'timefreq';
                            'plot'     'string' { 'no','alpha','sigout','sigoutm','sigoutp','sigoutm2' }  'sigoutm'; 
                            'nofig'    'string' { 'on','off' } 'off';
                            'cbar'     'string' { 'on','off' } 'on';
                            'smooth'   'string' { 'on','off' } 'off';
                            'erpout'   'real'   []             [];
                            'tfopt'    'cell'   []             {};
                            'erpimopt' 'cell'   []             {} });
if isstr(g), error(g); end;

fprintf('Generating %d frequencies in log scale (ignore message on linear scale)\n', g.freqs(2));
g.freqs = logscale(g.freqs(1), g.freqs(3), g.freqs(2));
    
frames = length(timevect);
if size(data,1) == 1
    data = reshape(data, frames, size(data,2)*size(data,3)/frames);
end;

% trim sortvar values
% -------------------
[sortvar sortorder] = sort(sortvar);
data = data(:, sortorder);
len = length(sortvar);
lowindex  = round(len*g.trim(1)/100)+1;
highindex = round(len*g.trim(2)/100);
sortvar = sortvar(lowindex:highindex);
data    = data(:, lowindex:highindex);
if lowindex ~=1 | highindex ~= length(sortorder)
    fprintf('Actual percentiles %1.2f-%1.2f (indices 1-%d -> %d-%d): event vals min %3.2f; max %3.2f\n', ...
             100*(lowindex-1)/len, 100*highindex/len, len, lowindex, highindex, min(sortvar), max(sortvar));
end;

% assign subject number for each trial
% ------------------------------------
if ~isempty(g.triallimit)
    alltrials = zeros(1,len);
    for index = 1:length(g.triallimit)-1
        alltrials([g.triallimit(index):g.triallimit(index+1)-1]) = index; 
    end;
    alltrials = alltrials(sortorder);
    alltrials = alltrials(lowindex:highindex);
end;

%figure; hist(sortvar)

% constraining time values depending on data alignment
% ----------------------------------------------------
srate = 1000*(length(timevect)-1)/(timevect(end)-timevect(1));
if ~isempty(g.align)

    if iscell(g.align)
        fprintf('Realigned sortvar plotted at %g ms.\n',g.align{1});
        shifts = round((g.align{1}-g.align{2})*srate/1000); % shifts can be positive or negative
        g.align = g.align{1};
    else
        if isinf(g.align)
            g.align = median(sortvar);
        end
        fprintf('Realigned sortvar plotted at %g ms.\n',g.align);
        
        % compute shifts
        % --------------
        shifts = round((g.align-sortvar)*srate/1000); % shifts can be positive or negative
    end;
    
    %figure; hist(shifts)
    minshift = min(shifts); % negative
    maxshift = max(shifts); % positive
    if minshift > 0, error('minshift has to be negative'); end;
    if maxshift < 0, error('maxshift has to be positive'); end;
    
    % realign data for all trials
    % ---------------------------
    aligndata=zeros(size(data))+NaN; % begin with matrix of zeros()
    for t=1:size(data,2), %%%%%%%%% foreach trial %%%%%%%%%
        shft = shifts(t);
        if shft>0, % shift right
            aligndata(shft+1:frames,t)=data(1:frames-shft,t);
        elseif shft < 0 % shift left
            aligndata(1:frames+shft,t)=data(1-shft:frames,t);
        else % shft == 0
            aligndata(:,t) = data(:,t);
        end 
    end % end trial    
    data = aligndata(maxshift+1:frames+minshift,:);
    if any(any(isnan(data))), error('NaNs remaining after data alignment'); end;
    timevect = timevect(maxshift+1:frames+minshift);
    
    % take the time vector subset
    % ---------------------------
    if isempty(timevect), error('Shift too big, empty time vector'); 
    else fprintf('Time vector truncated for data alignment between %1.0f and %1.0f\n', ...
                 min(timevect), max(timevect)); 
    end;        
end;

% regress out the ERP from the data (4 times to have residual ERP close to 0)
% ---------------------------------------------------------------------------
if ~isempty(g.erpout) 
    %data = erpregout(data);
    %data = erpregout(data);
    %data = erpregout(data);
    disp('Regressing out ERP');
    erpbeg = mean(data,2);
    if ~isempty(g.triallimit)
        for index = 1:length(g.triallimit)-1
            trials = find(alltrials == index); 
            %[data(:,trials) erp factors]= erpregout(data(:,trials), [timevect(1) timevect(end)], [300 400]);
            erpstart = mean(data(:,trials),2);
            data(:,trials) = erpregout(data(:,trials));
            data(:,trials) = erpregout(data(:,trials));
            data(:,trials) = erpregout(data(:,trials));
            data(:,trials) = erpregout(data(:,trials));
            erpend = mean(data(:,trials),2);
            fprintf([ '***********************************************\n' ...
                      'Ratio in ERP after regression (compare to before) is %3.2f\n' ...
                      '***********************************************\n' ], mean(erpend./erpstart));
        end;
    end;
    erpend = mean(data,2);
    fprintf([ '***********************************************\n' ...
              'Ratio in grand ERP after regression (compare to before) is %3.2f\n' ...
              '***********************************************\n' ], mean(erpend./erpbeg));
    
    if 0
    %trials = find(alltrials == 1); 
    data2 = data;
    dasf

    [data(:,trials) erp factors]= erpregout(data(:,trials), [timevect(1) timevect(end)], [300 400]);

    figure; subplot(1,2,1); erpimage(data, sortvar, timevect, '', 300, 500, 'erp');
    
    figure; 
    subplot(1,2,1); erpimage(data2, sortvar, timevect, '', 300, 500, 'erp');
    subplot(1,2,2); erpimage(data , sortvar, timevect, '', 300, 500, 'erp');
    
    figure; 
    subplot(1,2,1); erpimage(data2(:,trials), sortvar(trials), timevect, '', 0, 0, 'erp');
    subplot(1,2,2); erpimage(data (:,trials), sortvar(trials), timevect, '', 0, 0, 'erp');

    
    figure; 
    for index = 1:length(g.triallimit)-1
        subplot(3,5,index);
        trials = find(alltrials == index); 
        erpimage(data(:,trials), sortvar(trials), timevect, '', 50, 100, 'erp');
    end;
    
    disp('Regressing out ERP');
    if ~isempty(g.triallimit)
        for index = 1:length(g.triallimit)-1
            trials = find(alltrials == index); 
            [data(:,trials) erp factors]= erpregout(data(:,trials), [timevect(1) timevect(end)], [300 400]);
            data(:,trials) = erpregout(data(:,trials));
            data(:,trials) = erpregout(data(:,trials));
            data(:,trials) = erpregout(data(:,trials));
            data(:,trials) = erpregout(data(:,trials));
            data(:,trials) = erpregout(data(:,trials));
            data(:,trials) = erpregout(data(:,trials));
            data(:,trials) = erpregout(data(:,trials));
            %data(:,trials) = erpregout(data(:,trials), [timevect(1) timevect(end)], g.erpout);
            %data(:,trials) = erpregout(data(:,trials), [timevect(1) timevect(end)], g.erpout);
            %data(:,trials) = erpregout(data(:,trials), [timevect(1) timevect(end)], g.erpout);
            %data(:,trials) = erpregout(data(:,trials), [timevect(1) timevect(end)], g.erpout);
            %data(:,trials) = erpregout(data(:,trials), [timevect(1) timevect(end)], g.erpout);
            %data(:,trials) = erpregout(data(:,trials), [timevect(1) timevect(end)], g.erpout);
            %data(:,trials) = erpregout(data(:,trials), [timevect(1) timevect(end)], g.erpout);
            %data(:,trials) = erpregout(data(:,trials), [timevect(1) timevect(end)], g.erpout);
        end;
    else
        %data = erpregout(data, [timevect(1) timevect(end)], g.erpout);
    end;
    end;
end;

% generate time points
% --------------------
g.times = gettimes(timevect, g.times);
data = reshape(data, 1, size(data,2)*size(data,1));

% time frequency decomposition
% ----------------------------
if strcmpi(g.method, 'timefreq') & isempty(g.plotvals)
    data = reshape(data, length(data)/length(sortvar), length(sortvar));
    [tf, g.freqs, g.times] = timefreq(data, srate, 'freqs', g.freqs, 'timesout', g.times, ...
                                      'tlimits', [min(timevect) max(timevect)], 'wavelet', 3);
    outvar   = sortvar;
end;

% compute correlation
% -------------------
if strcmpi(g.mode, 'phase') & isempty(g.plotvals)    
    for freq = 1:length(g.freqs)
        fprintf('Computing correlation with phase %3.2f Hz ----------------------\n', g.freqs(freq));
        for time = 1:length(g.times)
            if strcmpi(g.method, 'erpimage')
                [outdata,outvar,outtrials,limits,axhndls,erp, ...
                 amps,cohers,cohsig,ampsig,outamps,phsangls,phsamp,sortidx] =erpimage(data,sortvar,timevect, ...
                        '', 0,0,g.erpimopt{:}, 'phasesort', [g.times(time) 0 g.freqs(freq)], 'noshow', 'yes');
            else
                phsangls = angle( squeeze(tf(freq, time, :))' );
            end;
            
            % computing ITCs
            [ ITC(freq, time) alpha(freq, time) ] = ...
                    bootcircle(outvar/mean(outvar), exp(j*phsangls), 'naccu', 250);

            % accumulate 200 values, fitted with a normal distribution
            % --------------------------------------------------------
            %cmplx = outvar .* exp(j*phsangls)/mean(outvar);
            %ITC(freq, time) = mean(cmplx);
            %alpha(freq,time) = bootstat(outvar/mean(outvar), exp(j*phsangls), 'res = mean(arg1 .* arg2)', ...
            %                         'naccu', 100, 'vals', abs(ITC(freq, time)), 'distfit', 'norm');
        end;
    end;
    sigout = ITC;
elseif isempty(g.plotvals)
    for freq = 1:length(g.freqs)
        fprintf('Computing correlation with amplitude %3.2f Hz ----------------------\n', g.freqs(freq));
        for time = 1:length(g.times)
            if strcmpi(g.method, 'erpimage')
                [outdata,outvar,outtrials,limits,axhndls,erp, ...
                 amps,cohers,cohsig,ampsig,outamps,phsangls,phsamp,sortidx] =erpimage(data,sortvar,timevect, ...
                             '', 0,0,g.erpimopt{:}, 'ampsort', [g.times(time) 0 g.freqs(freq)], 'noshow', 'yes');
            else
                phsamp = abs(squeeze(tf(freq, time, :)))';
            end;
        
            % computing ITCs
            [ypred alpha(freq, time) Rsq slope(freq, time)] = myregress(outvar, 20*log10(phsamp));
        end;
    end;    
    sigout = slope;
else 
    g.times = g.plotvals{1};
    g.freqs = g.plotvals{2};
    alpha   = g.plotvals{3};
    sigout  = g.plotvals{4};
end;

% plot correlation
% ----------------
if ~strcmp('plot', 'no')
    if ~isreal(sigout), 
        if strcmpi(g.plot, 'sigoutp')
            sigoutplot = angle(sigout);
        else
            sigoutplot = abs(sigout);
        end;
    else sigoutplot = sigout; 
    end;
    sigouttmp = sigoutplot;
    if strcmpi(g.smooth, 'on')
        tmpfilt = gauss2d(3,3,.3,.3);
        tmpfilt = tmpfilt/sum(sum(tmpfilt));
        alpha = conv2(alpha, tmpfilt, 'same');
    end;
    
    % mask signal out
    if g.pmask > 0.5
        indices = find( alpha(:) < g.pmask);
        sigouttmp = sigoutplot;
        sigouttmp(indices) = 0;
    elseif g.pmask < 0 % both sides, i.e. 0.01
        sigouttmp = sigoutplot;
        indices = intersect(find( alpha(:) > -g.pmask), find( alpha(:) < 1+g.pmask));
        sigouttmp(indices) = 0;
    else
        sigouttmp = sigoutplot;
        indices = find( alpha(:) > g.pmask);
        sigouttmp(indices) = 0;
    end;
    
    if strcmpi(g.nofig, 'off'), figure; end;
    switch g.plot
     case 'alpha',   limits = plotfig(g.times, g.freqs, -log10(alpha), g);
     case 'sigout',  limits = plotfigsim(g.times, g.freqs, sigoutplot, g);
     case { 'sigoutm' 'sigoutp' }, limits = plotfigsim(g.times, g.freqs, sigouttmp, g);
     case 'sigoutm2', limits = subplot(1,2,1); plotfigsim(g.times, g.freqs, sigoutplot, g);
                      limits = subplot(1,2,2); plotfigsim(g.times, g.freqs, sigouttmp, g);
    end;
end;

time = g.times;
freq = g.freqs;

return;

% formula for the normal distribution
% -----------------------------------
function y = norm(mu, sigma, x);

y = 1/sqrt(2) * exp( -(x-mu).^2/(sigma*sigma*2) ) / (sqrt(pi)*sigma);

% get time points
% ---------------
function val = gettimes(timevect, timevar);
    if length(timevar) == 2

        if timevar(1) > 0
            % generate linearly space vector
            % ------------------------------
            npoints = timevar(1);
            trim    = timevar(2);
            if length(timevect)-2*round(trim/100*length(timevect)) < npoints
                npoints = length(timevect)-round(trim/100*length(timevect));
            end;
            fprintf('Generating %d times points trimmed by %1.1f percent\n', npoints, trim);
            timer = max(timevect) - min(timevect);
            maxt = max(timevect)-timer*trim/100;
            mint = min(timevect)+timer*trim/100;
            val = linspace(mint,maxt, npoints);
        else            
            % subsample data
            % --------------
            nsub    = -timevar(1);
            trim    = timevar(2);
            len     = length(timevect);
            trimtimevect = timevect(round(trim/100*len)+1:len-round(trim/100*len));
            fprintf('Subsampling by %d and triming data by %1.1f percent (%d points)\n', nsub, trim, round(trim/100*len));
            val = trimtimevect(1:nsub:end);
        end;
    else
        val = timevar;
    end;        
    
    % find closet points in data
    oldval = val;
    for index = 1:length(val)
        [dum ind] = min(abs(timevect-val(index)));
        val(index) = timevect(ind);
    end;
    if length(val) < length(unique(val))
        disp('Warning: duplicate times, reduce the number of output times');
    end;
    if all(oldval == val)
        disp('Debug msg: Time value unchanged by finding closest in data');
    end;    
    
% get log scale (for frequency)
% ------------- 
function val = logscale(a,b,n);
    %val = [5 7 9]; return;
    val = linspace(log(a), log(b), n);
    val = exp(val);
    
% plot figure
% -----------
function limits = plotfig(times, freqs, vals, g)
    
    imagesc(times, [1:size(vals,1)], vals);

    ticks = linspace(1, size(vals,1), length(freqs));
    ticks = ticks(1:4:end);
    set(gca, 'ytick', ticks);
    set(gca, 'yticklabel', num2str(freqs(1:4:end)', 3))

    xlabel('Time (ms)'); ylabel('Freq (Hz)');
    if ~isempty(g.limits), caxis(g.limits); end;
    limits = caxis;
    if ~isempty(g.vert)
        for vert = g.vert(:)'
            hold on; plot([vert vert], [0.001 500], 'k', 'linewidth', 2);
        end;
    end;
    if strcmpi(g.cbar,'on')
        cbar;
    end;

% plot figure with symetrical phase
% ---------------------------------
function limits = plotfigsim(times, freqs, vals, g)
    
    imagesc(times, [1:size(vals,1)], vals);

    ticks = linspace(1, size(vals,1), length(freqs));
    ticks = ticks(1:4:end);
    set(gca, 'ytick', ticks);
    set(gca, 'yticklabel', num2str(freqs(1:4:end)', 3))

    if ~isempty(g.limits) 
        caxis(g.limits);
        limits = g.limits;
    else 
        clim = max(abs(caxis));
        limits = [-clim clim];
        caxis(limits);
    end;
    xlabel('Time (ms)'); ylabel('Freq (Hz)');
    if ~isempty(g.vert)
        for vert = g.vert(:)'
            hold on; plot([vert vert], [0.01 500], 'k', 'linewidth', 2);
        end;
    end;
    if strcmpi(g.cbar,'on')
        cbar;
    end;
    return;


% -----------------------------------------
% plot time-freq point (when clicking on it
% -----------------------------------------
function plotpoint(data, sortvar, timevect, freq, timepnts);

figure; 
subplot(2,2,1);
erpimage(act_all(:,:),sortvar_all,timepnts, ...
   '', 300,10,'erp', 'erpstd', 'caxis',[-1.0 1.0], 'srate',256,'align',352, 'yerplabel', '', erpimopt{:}, ...
                                                   'phasesort', [500 0 5]); % aligned to median rt

% plot in polar coordinates phase versus RT
% -----------------------------------------
phaseang2 = [phsangls phsangls-2*pi]; phaseang2 = movav(phaseang2,[], 100); 
outvar2   = [outvar   outvar];   outvar2   = movav(outvar2,[], 100);
phaseang2 = phaseang2(length(phsangls)/2-50:length(phsangls)+length(phsangls)/2-50);
outvar2   = outvar2  (length(phsangls)/2-50:length(phsangls)+length(phsangls)/2-50);
subplot(2,2,3);
polar(phsangls, outvar, '.');
hold on; polar(phaseang2, outvar2, 'r');

% computing ITC
% -------------
cmplx = outvar .* exp(j*phsangls);
ITCval = mean(cmplx);
angle(ITCval)/pi*180;
abs(ITCval)/mean(outvar);
text(-1300,-1400,sprintf('ITC value: amplitude %3.4f, angle %3.1f\n', abs(ITCval)/mean(outvar), angle(ITCval)/pi*180));

% accumulate 200 values
% ---------------------
alpha = 0.05;
if exist('accarray') ~= 1, accarray = NaN; end;
[sigval accarray] = bootstat(outvar/mean(outvar), phsangls, 'res = mean(arg1 .* exp(arg2*complex(0,1)))', ...
                      'accarray', accarray, 'bootside', 'upper', 'alpha', alpha);
text(-1300,-1600,sprintf('Threshold for significance at %d percent(naccu=200): %3.4f\n', alpha*100, sigval));
title(sprintf('Clust %d corr. theta phase at 500 ms and RT', clust));

% amplitude sorting
% -----------------
figure;
[outdata,outvar,outtrials,limits,axhndls,erp, ...
          amps,cohers,cohsig,ampsig,outamps,phsangls,phsamp,sortidx]=erpimage(act_all(:,:),sortvar_all,timepnts, ...
   '', 0,0,'erp', 'erpstd', 'caxis', 0.5, 'cbar', ...
   'srate',256,'align',352, 'yerplabel', '', erpimopt{:}, 'ampsort', [500 0 5]); % aligned to median rt
close;
subplot(2,2,2);
erpimage(act_all(:,:),sortvar_all,timevect, ...
   '', 300,10,'erp', 'erpstd', 'caxis', 0.5, 'cbar', ...
   'srate',256,'align',352, 'yerplabel', '', erpimopt{:}, 'ampsort', [500 0 5]); % aligned to median rt

% compute correlation
% -------------------
[ypred alpha Rsq] = myregress(outvar, phsamp);
subplot(2,2,4);
plot(outvar, phsamp, '.');
hold on;
plot(outvar, ypred, 'r');
xlabel('Reaction time');
ylabel('Amplitude');
title(sprintf('Theta amp. at 500 ms vs RT (p=%1.8f, R^2=%3.4f)', alpha, Rsq));
set(gcf, 'position', [336 485 730 540], 'paperpositionmode', 'auto');

setfont(gcf, 'fontsize', 10);
eval(['print -djpeg clust' int2str(clust) 'corrthetart.jpg']);
    
    
% set up for command line call
% copy text from plotcorrthetaart
%timevect = times;
sortvar = sortvar_all;
data = act_all;
time = 1
freq = 1
g.times = 0;
g.freqs = 5;
g.erpimopt = {};

