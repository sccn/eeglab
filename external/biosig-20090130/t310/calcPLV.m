function r = calcPLV(s, h, startend, f, trls, class, rmartif, alpha)
% Calculates PLV maps.
%
% This function calculates phase locking value (PLV) maps by using Gabor
% wavelets to estimate the synchronization between two electrode sites in
% specific frequency bands.
%
% Usage:
%   r = calcPLV(s, h, startend, f, trls, class, rmartif, alpha);
%
% Input parameters:
%   s        ... Input signal (as obtained by sload), must have exactly two
%                channels (columns)
%   h        ... Header structure (as obtained by sload)
%   startend ... Start and end point within trial (s)
%   f        ... Frequency range (Hz)
%
% Optional input parameters:
%   trls    ... Trials
%   class   ... Class
%   rmartif ... Remove artifacts
%   alpha   ... Significance level
%
% Output parameter:
%   r ... Structure containing the results

% Copyright by Clemens Brunner
% $Revision: 1.1 $ $Date: 2009-01-30 06:04:44 $
% E-Mail: clemens.brunner@tugraz.at

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

if (nargin < 1)
    error('No input signal specified.');
end;
if (nargin < 2)
    error('No header structure specified.');
end;
if (nargin < 3)
    error('Start and end point of one trial not specified.');
end;
if (nargin < 4)
    error('Frequency range not specified.');
end;
if (nargin < 5)
    trls = [];
end;
if (nargin < 6)
    class = [];  % All classes
end;
if (nargin < 7)
    rmartif = false;
end;
if (nargin < 8)
    alpha = 0.05;
end;

if (isempty(class))
    class = unique(h.Classlabel);
end;

% Trigger data
fs = h.SampleRate;

if (rmartif)
    [st, temp] = trigg(s, h.TRIG(ismember(h.Classlabel, class) & ...
                       h.ArtifactSelection == 0), startend(1)*fs, ...
                       startend(2)*fs-1);
else
    [st, temp] = trigg(s, h.TRIG(ismember(h.Classlabel, class)), ...
                       startend(1)*fs, startend(2)*fs-1);
end;

% Select specified trials
if ~isempty(trls)
    st = reshape(st, temp);
    st = st(:, :, trls);
    st = reshape(st, temp(1), temp(2)*length(trls));
end;

triallen = (startend(2) - startend(1)) * fs;
trials = length(st) / triallen;

st = st';

N = 1;  % Number of surrogate data series
alpha = alpha/(triallen * length(f));  % Confidence level with Bonferroni correction
plv = zeros(triallen, length(f));
sig = zeros(length(f), triallen);

% Calculates N PLVs between original signal 1 and shuffled signal 2 for all time
% instants and all frequencies
for m = 1:length(f)
    PLSn = zeros(triallen, N);
    disp(['f = ' num2str(f(m)) 'Hz']);
    for n = 1:N
        disp(['  ' num2str(100*(n-1)/N) '%']);
        temp = gaborPLV(st(:,1), shuffleTrials(st(:,2), trials), fs, trials, triallen, f(m));
        PLSn(:,n) = temp;  % Each column corresponds to the PLV of a particular surrogate data set 
    end;
    
    plv(:,m) = gaborPLV(st(:,1), st(:,2), fs, trials, triallen, f(m));  % PLV of original signals
    temp = repmat(plv(:,m)',N,1);  % Matrix is used for comparison only
    sig(m,:) = sum(temp<PLSn')./N;  % How many values of PLSn > original PLV (in %)?
    sig(m,:) = sig(m,:)<=alpha;  % Significant, if <= alpha
end;

sig(sig==0) = NaN;
sig=sig';

r{1}.plv = plv;
r{1}.sig = sig;
r{1}.fs = fs;
r{1}.f = f;
r{1}.datatype = 'PLV';