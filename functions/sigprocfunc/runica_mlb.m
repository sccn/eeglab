% runica() - Perform Independent Component Analysis (ICA) decomposition
%            of input data using the logistic infomax ICA algorithm of
%            Bell & Sejnowski (1995) with the natural gradient feature
%            of Amari, Cichocki & Yang, or optionally the extended-ICA
%            algorithm of Lee, Girolami & Sejnowski, with optional PCA
%            dimension reduction. Annealing based on weight changes is
%            used to automate the separation process.
% Usage:
%         >> [weights,sphere] = runica(data); % train using defaults
%    else
%         >> [weights,sphere,compvars,bias,signs,lrates,activations] ...
%                             = runica(data,'Key1',Value1',...);
% Input:
%    data     = input data (chans,frames*epochs).
%               Note that if data consists of multiple discontinuous epochs,
%               each epoch should be separately baseline-zero'd using
%                  >> data = rmbase(data,frames,basevector);
%
% Optional keywords [argument]:
% 'extended'  = [N] perform tanh() "extended-ICA" with sign estimation
%               N training blocks. If N > 0, automatically estimate the
%               number of sub-Gaussian sources. If N < 0, fix number of
%               sub-Gaussian comps to -N [faster than N>0] (default|0 -> off)
% 'pca'       = [N] decompose a principal component     (default -> 0=off)
%               subspace of the data. Value is the number of PCs to retain.
% 'ncomps'    = [N] number of ICA components to compute (default -> chans or 'pca' arg)
%               using rectangular ICA decomposition
% 'sphering'  = ['on'/'off'] flag sphering of data      (default -> 'on')
% 'weights'   = [W] initial weight matrix               (default -> eye())
%                            (Note: if 'sphering' 'off', default -> spher())
% 'lrate'     = [rate] initial ICA learning rate (<< 1) (default -> heuristic)
% 'block'     = [N] ICA block size (<< datalength)      (default -> heuristic)
% 'anneal'    = annealing constant (0,1] (defaults -> 0.90, or 0.98, extended)
%                         controls speed of convergence
% 'annealdeg' = [N] degrees weight change for annealing (default -> 70)
% 'stop'      = [f] stop training when weight-change < this (default -> 1e-6
%               if less than 33 channel and 1E-7 otherwise)
% 'maxsteps'  = [N] max number of ICA training steps    (default -> 512)
% 'bias'      = ['on'/'off'] perform bias adjustment    (default -> 'on')
% 'momentum'  = [0<f<1] training momentum               (default -> 0)
% 'specgram'  = [srate loHz hiHz frames winframes] decompose a complex time/frequency
%               transform of the data (Note: winframes must divide frames)
%                            (defaults [srate 0 srate/2 size(data,2) size(data,2)])
% 'posact'    = make all component activations net-positive(default 'on'}
% 'verbose'   = give ascii messages ('on'/'off')        (default -> 'on')
%
% Outputs:    [Note: RO means output in reverse order of projected mean variance
%                    unless starting weight matrix passed ('weights' above)]
% weights     = ICA weight matrix (comps,chans)      [RO]
% sphere      = data sphering matrix (chans,chans) = spher(data)
%               Note that unmixing_matrix = weights*sphere {if sphering off -> eye(chans)}
% compvars    = back-projected component variances   [RO]
% bias        = vector of final (ncomps) online bias [RO]    (default = zeros())
% signs       = extended-ICA signs for components    [RO]    (default = ones())
%                   [ -1 = sub-Gaussian; 1 = super-Gaussian]
% lrates      = vector of learning rates used at each training step [RO]
% activations = activation time courses of the output components (ncomps,frames*epochs)
%
% Authors: Scott Makeig with contributions from Tony Bell, Te-Won Lee,
% Tzyy-Ping Jung, Sigurd Enghoff, Michael Zibulevsky, Delorme Arnaud,
% CNL/The Salk Institute, La Jolla, 1996-

% Uses: posact()

% Reference (please cite):
%
% Makeig, S., Bell, A.J., Jung, T-P and Sejnowski, T.J.,
% "Independent component analysis of electroencephalographic data,"
% In: D. Touretzky, M. Mozer and M. Hasselmo (Eds). Advances in Neural
% Information Processing Systems 8:145-151, MIT Press, Cambridge, MA (1996).
%
% Toolbox Citation:
%
% Makeig, Scott et al. "EEGLAB: ICA Toolbox for Psychophysiological Research".
% WWW Site, Swartz Center for Computational Neuroscience, Institute of Neural
% Computation, University of San Diego California
% <www.sccn.ucsd.edu/eeglab/>, 2000. [World Wide Web Publication].
%
% For more information:
% http://www.sccn.ucsd.edu/eeglab/icafaq.html - FAQ on ICA/EEG
% http://www.sccn.ucsd.edu/eeglab/icabib.html - mss. on ICA & biosignals
% http://www.cnl.salk.edu/~tony/ica.html - math. mss. on ICA

% Copyright (C) 1996 Scott Makeig et al, SCCN/INC/UCSD, scott@sccn.ucsd.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

%%%%%%%%%%%%%%%%%%%%%%%%%%% Edit history %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  runica()  - by Scott Makeig with contributions from Tony Bell, Te-Won Lee
%              Tzyy-Ping Jung, Sigurd Enghoff, Michael Zibulevsky et al.
%                            CNL / Salk Institute 1996-00
%  04-30-96 built from icatest.m and ~jung/.../wtwpwica.m -sm
%  07-28-97 new runica(), adds bias (default on), momentum (default off),
%           extended-ICA (Lee & Sejnowski, 1997), cumulative angledelta
%           (until lrate drops), keywords, signcount for speeding extended-ICA
%  10-07-97 put acos() outside verbose loop; verbose 'off' wasn't stopping -sm
%  11-11-97 adjusted help msg -sm
%  11-30-97 return eye(chans) if sphering 'off' or 'none' (undocumented option) -sm
%  02-27-98 use pinv() instead of inv() to rank order comps if ncomps < chans -sm
%  04-28-98 added 'posact' and 'pca' flags  -sm
%  07-16-98 reduced length of randperm() for kurtosis subset calc. -se & sm
%  07-19-98 fixed typo in weights def. above -tl & sm
%  12-21-99 added 'specgram' option suggested by Michael Zibulevsky, UNM -sm
%  12-22-99 fixed rand() sizing inefficiency on suggestion of Mike Spratling, UK -sm
%  01-11-00 fixed rand() sizing bug on suggestion of Jack Foucher, Strasbourg -sm
%  12-18-00 test for existence of Sig Proc Tlbx function 'specgram'; improve
%           'specgram' option arguments -sm
%  01-25-02 reformated help & license -ad
%  01-25-02 lowered default lrate and block -ad
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [weights,sphere,meanvar,bias,signs,lrates,activations,y,loglik] = runica(data,p1,v1,p2,v2,p3,v3,p4,v4,p5,v5,p6,v6,p7,v7,p8,v8,p9,v9,p10,v10,p11,v11,p12,v12,p13,v13,p14,v14)

if nargin < 1
    help runica
    return
end

[chans frames] = size(data); % determine the data size
urchans = chans;  % remember original data channels
datalength = frames;
if chans<2
    fprintf('\nrunica() - data size (%d,%d) too small.\n\n', chans,frames);
    return
end
%
%%%%%%%%%%%%%%%%%%%%%% Declare defaults used below %%%%%%%%%%%%%%%%%%%%%%%%
%
MAX_WEIGHT           = 1e8;       % guess that weights larger than this have blown up
DEFAULT_STOP         = 1e-6;  % stop training if weight changes below this
DEFAULT_ANNEALDEG    = 60;        % when angle change reaches this value,
DEFAULT_ANNEALSTEP   = 0.98;      %     anneal by multiplying lrate by this
DEFAULT_EXTANNEAL    = 0.98;      %     or this if extended-ICA
DEFAULT_MAXSTEPS     = 500;       % stop training after this many steps
DEFAULT_MOMENTUM     = 0.0;       % default momentum weight

DEFAULT_BLOWUP       = 1000000000.0;   % = learning rate has 'blown up'
DEFAULT_BLOWUP_FAC   = 0.8;       % when lrate 'blows up,' anneal by this fac
DEFAULT_RESTART_FAC  = 0.9;       % if weights blowup, restart with lrate
% lower by this factor
MIN_LRATE            = 0.0001;  % if weight blowups make lrate < this, quit
MAX_LRATE            = 0.1;       % guard against uselessly high learning rate
DEFAULT_LRATE        = 0.00065/log(chans);
% heuristic default - may need adjustment
%   for large or tiny data sets!
% DEFAULT_BLOCK        = floor(sqrt(frames/4));  % heuristic default
DEFAULT_BLOCK          = ceil(min(5*log(frames),0.3*frames)); % heuristic
% - may need adjustment!
% Extended-ICA option:
DEFAULT_EXTENDED     = 1;         % default off
DEFAULT_EXTBLOCKS    = 1;         % number of blocks per kurtosis calculation
DEFAULT_NSUB         = 0;         % initial default number of assumed sub-Gaussians
% for extended-ICA
DEFAULT_EXTMOMENTUM  = 0.5;       % momentum term for computing extended-ICA kurtosis
MAX_KURTSIZE         = 6000;      % max points to use in kurtosis calculation
MIN_KURTSIZE         = 2000;      % minimum good kurtosis size (flag warning)
SIGNCOUNT_THRESHOLD  = 25;        % raise extblocks when sign vector unchanged
% after this many steps
SIGNCOUNT_STEP       = 2;         % extblocks increment factor

DEFAULT_SPHEREFLAG   = 'on';      % use the sphere matrix as the default
%   starting weight matrix
DEFAULT_PCAFLAG      = 'off';     % don't use PCA reduction
DEFAULT_POSACTFLAG   = 'on';      % use posact()
DEFAULT_VERBOSE      = 1;         % write ascii info to calling screen
DEFAULT_BIASFLAG     = 1;         % default to using bias in the ICA update rule
%
%%%%%%%%%%%%%%%%%%%%%%% Set up keyword default values %%%%%%%%%%%%%%%%%%%%%%%%%
%
if nargout < 2,
    fprintf('runica() - needs at least two output arguments.\n');
    return
end
epochs = 1;							 % do not care how many epochs in data

pcaflag    = DEFAULT_PCAFLAG;
sphering   = DEFAULT_SPHEREFLAG;     % default flags
posactflag = DEFAULT_POSACTFLAG;
verbose    = DEFAULT_VERBOSE;

block      = DEFAULT_BLOCK;          % heuristic default - may need adjustment!
lrate      = DEFAULT_LRATE;
annealdeg  = DEFAULT_ANNEALDEG;
annealstep = 0;                      % defaults declared below
nochange   = NaN;
momentum   = DEFAULT_MOMENTUM;
maxsteps   = DEFAULT_MAXSTEPS;

weights    = 0;                      % defaults defined below
ncomps     = chans;
biasflag   = DEFAULT_BIASFLAG;

extended   = DEFAULT_EXTENDED;
extblocks  = DEFAULT_EXTBLOCKS;
kurtsize   = MAX_KURTSIZE;
signsbias  = 0.02;                   % bias towards super-Gaussian components
extmomentum= DEFAULT_EXTMOMENTUM;    % exp. average the kurtosis estimates
nsub       = DEFAULT_NSUB;
wts_blowup = 0;                      % flag =1 when weights too large
wts_passed = 0;                      % flag weights passed as argument
%
%%%%%%%%%% Collect keywords and values from argument list %%%%%%%%%%%%%%%
%
if (nargin> 1 && rem(nargin,2) == 0)
    fprintf('runica(): Even number of input arguments???')
    return
end
for i = 3:2:nargin % for each Keyword
    Keyword = eval(['p',int2str((i-3)/2 +1)]);
    Value = eval(['v',int2str((i-3)/2 +1)]);
    if ~ischar(Keyword)
        fprintf('runica(): keywords must be strings')
        return
    end
    Keyword = lower(Keyword); % convert upper or mixed case to lower

    if strcmp(Keyword,'weights') || strcmp(Keyword,'weight')
        if ischar(Value)
            fprintf(...
                'runica(): weights value must be a weight matrix or sphere')
            return
        else
            weights = Value;
            wts_passed =1;
        end
    elseif strcmp(Keyword,'ncomps')
        if ischar(Value)
            fprintf('runica(): ncomps value must be an integer')
            return
        end
        if ncomps < urchans && ncomps ~= Value
            fprintf('runica(): Use either PCA or ICA dimension reduction');
            return
        end
        ncomps = Value;
        if ~ncomps,
            ncomps = chans;
        end
    elseif strcmp(Keyword,'pca')
        if ncomps < urchans && ncomps ~= Value
            fprintf('runica(): Use either PCA or ICA dimension reduction');
            return
        end
        if ischar(Value)
            fprintf(...
                'runica(): pca value should be the number of principal components to retain')
            return
        end
        pcaflag = 'on';
        ncomps = Value;
        if ncomps > chans || ncomps < 1,
            fprintf('runica(): pca value must be in range [1,%d]\n',chans)
            return
        end
        chans = ncomps;
    elseif strcmp(Keyword,'posact')
        if ~ischar(Value)
            fprintf('runica(): posact value must be on or off')
            return
        else
            Value = lower(Value);
            if ~strcmp(Value,'on') && ~strcmp(Value,'off'),
                fprintf('runica(): posact value must be on or off')
                return
            end
            posactflag = Value;
        end
    elseif strcmp(Keyword,'lrate')
        if ischar(Value)
            fprintf('runica(): lrate value must be a number')
            return
        end
        lrate = Value;
        if lrate>MAX_LRATE || lrate <0,
            fprintf('runica(): lrate value is out of bounds');
            return
        end
        if ~lrate,
            lrate = DEFAULT_LRATE;
        end
    elseif strcmp(Keyword,'block') || strcmp(Keyword,'blocksize')
        if ischar(Value)
            fprintf('runica(): block size value must be a number')
            return
        end
        block = floor(Value);
        if ~block,
            block = DEFAULT_BLOCK;
        end
    elseif strcmp(Keyword,'stop') || strcmp(Keyword,'nochange') ...
            | strcmp(Keyword,'stopping')
        if ischar(Value)
            fprintf('runica(): stop wchange value must be a number')
            return
        end
        nochange = Value;
    elseif strcmp(Keyword,'maxsteps') || strcmp(Keyword,'steps')
        if ischar(Value)
            fprintf('runica(): maxsteps value must be an integer')
            return
        end
        maxsteps = Value;
        if ~maxsteps,
            maxsteps   = DEFAULT_MAXSTEPS;
        end
        if maxsteps < 0
            fprintf('runica(): maxsteps value (%d) must be a positive integer',maxsteps)
            return
        end
    elseif strcmp(Keyword,'anneal') || strcmp(Keyword,'annealstep')
        if ischar(Value)
            fprintf('runica(): anneal step value (%2.4f) must be a number (0,1)',Value)
            return
        end
        annealstep = Value;
        if annealstep <=0 || annealstep > 1,
            fprintf('runica(): anneal step value (%2.4f) must be (0,1]',annealstep)
            return
        end
    elseif strcmp(Keyword,'annealdeg') || strcmp(Keyword,'degrees')
        if ischar(Value)
            fprintf('runica(): annealdeg value must be a number')
            return
        end
        annealdeg = Value;
        if ~annealdeg,
            annealdeg = DEFAULT_ANNEALDEG;
        elseif annealdeg > 180 || annealdeg < 0
            fprintf('runica(): annealdeg (%3.1f) is out of bounds [0,180]',...
                annealdeg);
            return

        end
    elseif strcmp(Keyword,'momentum')
        if ischar(Value)
            fprintf('runica(): momentum value must be a number')
            return
        end
        momentum = Value;
        if momentum > 1.0 || momentum < 0
            fprintf('runica(): momentum value is out of bounds [0,1]')
            return
        end
    elseif strcmp(Keyword,'sphering') || strcmp(Keyword,'sphereing') ...
            | strcmp(Keyword,'sphere')
        if ~ischar(Value)
            fprintf('runica(): sphering value must be on, off, or none')
            return
        else
            Value = lower(Value);
            if ~strcmp(Value,'on') && ~strcmp(Value,'off') && ~strcmp(Value,'none'),
                fprintf('runica(): sphering value must be on or off')
                return
            end
            sphering = Value;
        end
    elseif strcmp(Keyword,'bias')
        if ~ischar(Value)
            fprintf('runica(): bias value must be on or off')
            return
        else
            Value = lower(Value);
            if strcmp(Value,'on')
                biasflag = 1;
            elseif strcmp(Value,'off'),
                biasflag = 0;
            else
                fprintf('runica(): bias value must be on or off')
                return
            end
        end
    elseif strcmp(Keyword,'specgram') || strcmp(Keyword,'spec')

        if ~exist('specgram') < 2 % if ~exist or defined workspace variable
            fprintf(...
                'runica(): MATLAB Sig. Proc. Toolbox function "specgram" not found.\n')
            return
        end
        if ischar(Value)
            fprintf('runica(): specgram argument must be a vector')
            return
        end
        srate = Value(1);
        if (srate < 0)
            fprintf('runica(): specgram srate (%4.1f) must be >=0',srate)
            return
        end
        if length(Value)>1
            loHz = Value(2);
            if (loHz < 0 || loHz > srate/2)
                fprintf('runica(): specgram loHz must be >=0 and <= srate/2 (%4.1f)',srate/2)
                return
            end
        else
            loHz = 0; % default
        end
        if length(Value)>2
            hiHz = Value(3);
            if (hiHz < loHz || hiHz > srate/2)
                fprintf('runica(): specgram hiHz must be >=loHz (%4.1f) and <= srate/2 (%4.1f)',loHz,srate/2)
                return
            end
        else
            hiHz = srate/2; % default
        end
        if length(Value)>3
            Hzframes = Value(5);
            if (Hzframes<0 || Hzframes > size(data,2))
                fprintf('runica(): specgram frames must be >=0 and <= data length (%d)',size(data,2))
                return
            end
        else
            Hzframes = size(data,2); % default
        end
        if length(Value)>4
            Hzwinlen = Value(4);
            if rem(Hzframes,Hzwinlen) % if winlen doesn't divide frames
                fprintf('runica(): specgram Hzinc must divide frames (%d)',Hzframes)
                return
            end
        else
            Hzwinlen = Hzframes; % default
        end
        Specgramflag = 1; % set flag to perform specgram()

    elseif strcmp(Keyword,'extended') || strcmp(Keyword,'extend')
        if ischar(Value)
            fprintf('runica(): extended value must be an integer (+/-)')
            return
        else
            extended = 1;      % turn on extended-ICA
            extblocks = fix(Value); % number of blocks per kurt() compute
            if extblocks < 0
                nsub = -1*fix(extblocks);  % fix this many sub-Gauss comps
            elseif ~extblocks,
                extended = 0;             % turn extended-ICA off
            elseif kurtsize>frames,   % length of kurtosis calculation
                kurtsize = frames;
                if kurtsize < MIN_KURTSIZE
                    fprintf(...
                        'runica() warning: kurtosis values inexact for << %d points.\n',...
                        MIN_KURTSIZE);
                end
            end
        end
    elseif strcmp(Keyword,'verbose')
        if ~ischar(Value)
            fprintf('runica(): verbose flag value must be on or off')
            return
        elseif strcmp(Value,'on'),
            verbose = 1;
        elseif strcmp(Value,'off'),
            verbose = 0;
        else
            fprintf('runica(): verbose flag value must be on or off')
            return
        end
    else
        fprintf('runica(): unknown flag')
        return
    end
end

%
%%%%%%%%%%%%%%%%%%%%%%%% Initialize weights, etc. %%%%%%%%%%%%%%%%%%%%%%%%
%
if ~annealstep,
    if ~extended,
        annealstep = DEFAULT_ANNEALSTEP;     % defaults defined above
    else
        annealstep = DEFAULT_EXTANNEAL;       % defaults defined above
    end
end % else use annealstep from commandline

if ~annealdeg,
    annealdeg  = DEFAULT_ANNEALDEG - momentum*90; % heuristic
    if annealdeg < 0,
        annealdeg = 0;
    end
end
if ncomps >  chans || ncomps < 1
    fprintf('runica(): number of components must be 1 to %d.\n',chans);
    return
end

if weights ~= 0,                    % initialize weights
    % starting weights are being passed to runica() from the commandline
    if verbose,
        fprintf('Using starting weight matrix named in argument list ...\n')
    end
    if  chans>ncomps && weights ~=0,
        [r,c]=size(weights);
        if r~=ncomps || c~=chans,
            fprintf(...
                'runica(): weight matrix must have %d rows, %d columns.\n', ...
                chans,ncomps);
            return;
        end
    end
end
%
%%%%%%%%%%%%%%%%%%%%% Check keyword values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if frames<chans,
    fprintf('runica(): data length (%d) < data channels (%d)!\n',frames,chans)
    return
elseif block < 1,
    fprintf('runica(): block size %d too small!\n',block)
    return
elseif block > frames,
    fprintf('runica(): block size exceeds data length!\n');
    return
elseif floor(epochs) ~= epochs,
    fprintf('runica(): data length is not a multiple of the epoch length!\n');
    return
elseif nsub > ncomps
    fprintf('runica(): there can be at most %d sub-Gaussian components!\n',ncomps);
    return
end

%
% adjust nochange if necessary
%
if isnan(nochange)
    if ncomps > 32
        nochange = 1E-7;
        nochangeupdated = 1; % for fprinting purposes
    else
        nochangeupdated = 1; % for fprinting purposes
        nochange = DEFAULT_STOP;
    end
else
    nochangeupdated = 0;
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Process the data %%%%%%%%%%%%%%%%%%%%%%%%%%
%
if verbose,
    fprintf( ...
        '\nInput data size [%d,%d] = %d channels, %d frames/n', ...
        chans,frames,chans,frames);
    if strcmp(pcaflag,'on')
        fprintf('After PCA dimension reduction,\n  finding ');
    else
        fprintf('Finding ');
    end
    if ~extended
        fprintf('%d ICA components using logistic ICA.\n',ncomps);
    else % if extended
        fprintf('%d ICA components using extended ICA.\n',ncomps);
        if extblocks > 0
            fprintf(...
                'Kurtosis will be calculated initially every %d blocks using %d data points.\n',...
                extblocks,     kurtsize);
        else
            fprintf(...
                'Kurtosis will not be calculated. Exactly %d sub-Gaussian components assumed.\n',...
                nsub);
        end
    end
    fprintf('Decomposing %d frames per ICA weight ((%d)^2 = %d weights, %d frames)\n',...
        floor(frames/ncomps.^2),ncomps.^2,frames);
    fprintf('Initial learning rate will be %g, block size %d.\n',...
        lrate,block);
    if momentum>0,
        fprintf('Momentum will be %g.\n',momentum);
    end
    fprintf( ...
        'Learning rate will be multiplied by %g whenever angledelta >= %g deg.\n', ...
        annealstep,annealdeg);

    if nochangeupdated
        fprintf('More than 32 channels: default stopping weight change 1E-7\n');
    end
    fprintf('Training will end when wchange < %g or after %d steps.\n', ...
        nochange,maxsteps);
    if biasflag,
        fprintf('Online bias adjustment will be used.\n');
    else
        fprintf('Online bias adjustment will not be used.\n');
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%% Remove overall row means %%%%%%%%%%%%%%%%%%%%%%%%
%
if verbose,
    fprintf('Removing mean of each channel ...\n');
end
data = data - mean(data')'*ones(1,frames);      % subtract row means

if verbose,
    fprintf('Final training data range: %g to %g\n', ...
        min(min(data)),max(max(data)));
end
%
%%%%%%%%%%%%%%%%%%% Perform PCA reduction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if strcmp(pcaflag,'on')
    fprintf('Reducing the data to %d principal dimensions...\n',ncomps);
    [eigenvectors,eigenvalues,data] = pcsquash(data,ncomps);
    % make data its projection onto the ncomps-dim principal subspace
end
%
%%%%%%%%%%%%%%%%%%% Perform specgram transformation %%%%%%%%%%%%%%%%%%%%%%%
%
if exist('Specgramflag') == 1
    % [P F T] = SPECGRAM(A,NFFT,Fs,WINDOW,NOVERLAP) % MATLAB Sig Proc Toolbox
    % Hzwinlen =  fix(srate/Hzinc); % CHANGED FROM THIS 12/18/00 -sm

    Hzfftlen = 2^(ceil(log(Hzwinlen)/log(2)));   % make FFT length next higher 2^k
    Hzoverlap = 0; % use sequential windows
    %
    % Get freqs and times from 1st channel analysis
    %
    [tmp,freqs,tms] = specgram(data(1,:),Hzfftlen,srate,Hzwinlen,Hzoverlap);

    fs = find(freqs>=loHz & freqs <= hiHz);
    if isempty(fs)
        fprintf('runica(): specified frequency range too narrow!\n');
        return
    end

    specdata = reshape(tmp(fs,:),1,length(fs)*size(tmp,2));
    specdata = [real(specdata) imag(specdata)];
    % fprintf('   size(fs) = %d,%d\n',size(fs,1),size(fs,2));
    % fprintf('   size(tmp) = %d,%d\n',size(tmp,1),size(tmp,2));
    %
    % Loop through remaining channels
    %
    for ch=2:chans
        [tmp] = specgram(data(ch,:),Hzwinlen,srate,Hzwinlen,Hzoverlap);
        tmp = reshape((tmp(fs,:)),1,length(fs)*size(tmp,2));
        specdata = [specdata;[real(tmp) imag(tmp)]]; % channels are rows
    end
    %
    % Print specgram confirmation and details
    %
    fprintf(...
        'Converted data to %d channels by %d=2*%dx%d points spectrogram data.\n',...
        chans,2*length(fs)*length(tms),length(fs),length(tms));
    if length(fs) > 1
        fprintf(...
            '   Low Hz %g, high Hz %g, Hz incr %g, window length %d\n',freqs(fs(1)),freqs(fs(end)),freqs(fs(2))-freqs(fs(1)),Hzwinlen);
    else
        fprintf(...
            '   Low Hz %g, high Hz %g, window length %d\n',freqs(fs(1)),freqs(fs(end)),Hzwinlen);
    end
    %
    % Replace data with specdata
    %
    data = specdata;
    datalength=size(data,2);
end
%
%%%%%%%%%%%%%%%%%%% Perform sphering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if strcmp(sphering,'on'), %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if verbose,
        fprintf('Computing the sphering matrix...\n');
    end
    sphere = 2.0*inv(sqrtm(cov(data'))); % find the "sphering" matrix = spher()
    if ~weights,
        if verbose,
            fprintf('Starting weights are the identity matrix ...\n');
        end
        weights = eye(ncomps,chans); % begin with the identity matrix
    else % weights given on commandline
        if verbose,
            fprintf('Using starting weights named on commandline ...\n');
        end
    end
    if verbose,
        fprintf('Sphering the data ...\n');
    end
    data = sphere*data;      % actually decorrelate the electrode signals

elseif strcmp(sphering,'off') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~weights
        if verbose,
            fprintf('Using the sphering matrix as the starting weight matrix ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
        sphere = 2.0*inv(sqrtm(cov(data'))); % find the "sphering" matrix = spher()
        weights = eye(ncomps,chans)*sphere; % begin with the identity matrix
        sphere = eye(chans);                 % return the identity matrix
    else % weights ~= 0
        if verbose,
            fprintf('Using starting weights named on commandline ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
        sphere = eye(chans);                 % return the identity matrix
    end
elseif strcmp(sphering,'none')
    sphere = eye(chans);                     % return the identity matrix
    if ~weights
        if verbose,
            fprintf('Starting weights are the identity matrix ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
        weights = eye(ncomps,chans); % begin with the identity matrix
    else % weights ~= 0
        if verbose,
            fprintf('Using starting weights named on commandline ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
    end
    sphere = eye(chans,chans);
    if verbose,
        fprintf('Returned variable "sphere" will be the identity matrix.\n');
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%% Initialize ICA training %%%%%%%%%%%%%%%%%%%%%%%%%
%
lastt=fix((datalength/block-1)*block+1);
BI=block*eye(ncomps,ncomps);
delta=zeros(1,chans*ncomps);
changes = [];
degconst = 180./pi;
startweights = weights;
prevweights = startweights;
oldweights = startweights;
prevwtchange = zeros(chans,ncomps);
oldwtchange = zeros(chans,ncomps);
lrates = zeros(1,maxsteps);
onesrow = ones(1,block);
bias = zeros(ncomps,1);
signs = ones(1,ncomps);    % initialize signs to nsub -1, rest +1
for k=1:nsub
    signs(k) = -1;
end
if extended && extblocks < 0 && verbose,
    fprintf('Fixed extended-ICA sign assignments:  ');
    for k=1:ncomps
        fprintf('%d ',signs(k));
    end; fprintf('\n');
end
signs = diag(signs); % make a diagonal matrix
oldsigns = zeros(size(signs));;
signcount = 0;              % counter for same-signs
signcounts = [];
urextblocks = extblocks;    % original value, for resets
old_kk = zeros(1,ncomps);   % for kurtosis momemtum
%
%%%%%%%% ICA training loop using the logistic sigmoid %%%%%%%%%%%%%%%%%%%
%
if verbose,
    fprintf('Beginning ICA training ...');
    if extended,
        fprintf(' first training step may be slow ...\n');
    else
        fprintf('\n');
    end
end
step=0;
laststep=0;
blockno = 1;  % running block counter for kurtosis interrupts
logstep = 1;  % iterator over log likelihood record
rand('state',sum(100*clock));  % set the random number generator state to
cost_step = 1; % record cost every cost_step iterations
% a position dependent on the system clock
while step < maxsteps, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    permute=randperm(datalength); % shuffle data order at each step
    loglik(logstep) = 0;
    for t=1:block:lastt, %%%%%%%%% ICA Training Block %%%%%%%%%%%%%%%%%%%
        pause(0);
        if ~isempty(get(0, 'currentfigure')) && strcmp(get(gcf, 'tag'), 'stop')
            close; error('USER ABORT');
        end
        if biasflag
            u=weights*data(:,permute(t:t+block-1)) + bias*onesrow;
        else
            u=weights*data(:,permute(t:t+block-1));
        end
        if ~extended
            %%%%%%%%%%%%%%%%%%% Logistic ICA weight update %%%%%%%%%%%%%%%%%%%
            y=1./(1+exp(-u));                                                %
            weights = weights + lrate*(BI+(1-2*y)*u')*weights;               %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else % Tanh extended-ICA weight update
            %%%%%%%%%%%%%%%%%%% Extended-ICA weight update %%%%%%%%%%%%%%%%%%%
            y=tanh(u);                                                       %
            weights = weights + lrate*(BI-signs*y*u'-u*u')*weights;          %
            
            %%%%%%%%% Calculate log likelihood given our model %%%%%%%%%
            if mod(step,cost_step) == 0
                subgcomp = find(diag(signs) == -1);
                supergcomp = find(diag(signs) == 1);
                loglik(logstep) = loglik(logstep) + sum(sum(log(exp(-(1/2)*(u(subgcomp,:)-1).^2) + exp(-(1/2)*(u(subgcomp,:)+1).^2))));
                loglik(logstep) = loglik(logstep) + sum(sum(-0.5*u(supergcomp,:).^2 - 2*log(cosh(u(supergcomp,:)))));
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        if biasflag
            if ~extended
                %%%%%%%%%%%%%%%%%%%%%%%% Logistic ICA bias %%%%%%%%%%%%%%%%%%%%%%%
                bias = bias + lrate*sum((1-2*y)')'; % for logistic nonlin. %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else % extended
                %%%%%%%%%%%%%%%%%%% Extended-ICA bias %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                bias = bias + lrate*sum((-2*y)')';  % for tanh() nonlin.   %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end

        if momentum > 0 %%%%%%%%% Add momentum %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            weights = weights + momentum*prevwtchange;
            prevwtchange = weights-prevweights;
            prevweights = weights;
        end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if max(max(abs(weights))) > MAX_WEIGHT
            wts_blowup = 1;
            change = nochange;
        end
        if extended && ~wts_blowup
            %
            %%%%%%%%%%% Extended-ICA kurtosis estimation %%%%%%%%%%%%%%%%%%%%%
            %
            if extblocks > 0 && rem(blockno,extblocks) == 0,
                % recompute signs vector using kurtosis
                if kurtsize < frames % 12-22-99 rand() size suggestion by M. Spratling
                    rp = fix(rand(1,kurtsize)*datalength);  % pick random subset
                    % Account for the possibility of a 0 generation by rand
                    ou = find(rp == 0);
                    while ~isempty(ou) % 1-11-00 suggestion by J. Foucher
                        rp(ou) = fix(rand(1,length(ou))*datalength);
                        ou = find(rp == 0);
                    end
                    partact=weights*data(:,rp(1:kurtsize));
                else                                        % for small data sets,
                    partact=weights*data;                   % use whole data
                end
                m2=mean(partact'.^2).^2;
                m4= mean(partact'.^4);
                kk= (m4./m2)-3.0;                           % kurtosis estimates
                if extmomentum
                    kk = extmomentum*old_kk + (1.0-extmomentum)*kk; % use momentum
                    old_kk = kk;
                end
                signs=diag(sign(kk+signsbias));             % pick component signs
                if signs == oldsigns,
                    signcount = signcount+1;
                else
                    signcount = 0;
                end
                oldsigns = signs;
                signcounts = [signcounts signcount];
                if signcount >= SIGNCOUNT_THRESHOLD,
                    extblocks = fix(extblocks * SIGNCOUNT_STEP);% make kurt() estimation
                    signcount = 0;                             % less frequent if sign
                end                                         % is not changing
            end % extblocks > 0 & . . .
        end % if extended %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        blockno = blockno + 1;
        if wts_blowup
            break
        end

    end % training block %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    loglik(logstep) = loglik(logstep) + log(abs(det(weights)));
    logstep = logstep + 1;
    
    if ~wts_blowup
        oldwtchange = weights-oldweights;
        step=step+1;
        %
        %%%%%%% Compute and print weight and update angle changes %%%%%%%%%
        %
        lrates(1,step) = lrate;
        angledelta=0.;
        delta=reshape(oldwtchange,1,chans*ncomps);
        change=delta*delta';
    end
    %
    %%%%%%%%%%%%%%%%%%%%%% Restart if weights blow up %%%%%%%%%%%%%%%%%%%%
    %
    if wts_blowup || isnan(change)|isinf(change),  % if weights blow up,
        fprintf('');
        step = 0;                          % start again
        change = nochange;
        wts_blowup = 0;                    % re-initialize variables
        blockno = 1;
        lrate = lrate*DEFAULT_RESTART_FAC; % with lower learning rate
        weights = startweights;            % and original weight matrix
        oldweights = startweights;
        change = nochange;
        oldwtchange = zeros(chans,ncomps);
        delta=zeros(1,chans*ncomps);
        olddelta = delta;
        extblocks = urextblocks;
        prevweights = startweights;
        prevwtchange = zeros(chans,ncomps);
        lrates = zeros(1,maxsteps);
        bias = zeros(ncomps,1);
        if extended
            signs = ones(1,ncomps);    % initialize signs to nsub -1, rest +1
            for k=1:nsub
                signs(k) = -1;
            end
            signs = diag(signs); % make a diagonal matrix
            oldsigns = zeros(size(signs));;
        end
        if lrate> MIN_LRATE
            r = rank(data);
            if r<ncomps
                fprintf('Data has rank %d. Cannot compute %d components.\n',...
                    r,ncomps);
                return
            else
                fprintf(...
                    'Lowering learning rate to %g and starting again.\n',lrate);
            end
        else
            fprintf( ...
                'runica(): QUITTING - weight matrix may not be invertible!\n');
            return;
        end
    else % if weights in bounds
        %
        %%%%%%%%%%%%% Print weight update information %%%%%%%%%%%%%%%%%%%%%%
        %
        if step> 2
            angledelta=acos((delta*olddelta')/sqrt(change*oldchange));
        end
        if verbose,
            places = -floor(log10(nochange));
            if step > 2,
                if ~extended,
                    ps = sprintf('step %d - lrate %5f, wchange %%%d.%df, angledelta %4.1f deg\n', ...
                        step,      lrate,        places+1,places,   degconst*angledelta);
                else
                    ps = sprintf('step %d - lrate %5f, wchange %%%d.%df, angledelta %4.1f deg, %d subgauss\n',...
                        step,      lrate,        degconst*angledelta,...
                        places+1,places,           (ncomps-sum(diag(signs)))/2);
                end
            elseif ~extended
                ps = sprintf('step %d - lrate %5f, wchange %%%d.%df\n',...
                    step,      lrate,       places+1,places  );
            else
                ps = sprintf('step %5d - lrate %5f, wchange %%%d.%df, %d subgauss\n',...
                    step,      lrate,        places+1,places, (ncomps-sum(diag(signs)))/2);
            end % step > 2
            fprintf('step %d - lrate %5f, wchange %8.8f, angledelta %5.1f deg, loglik %6.2f, nsub = %d\n', ...
                step,      lrate,     change, degconst*angledelta, loglik(step), sum(diag(signs)==-1));
            % fprintf(ps,change);  % <---- BUG !!!!
        end; % if verbose
        %
        %%%%%%%%%%%%%%%%%%%% Save current values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        changes = [changes change];
        oldweights = weights;
        %
        %%%%%%%%%%%%%%%%%%%% Anneal learning rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if degconst*angledelta > annealdeg,
            lrate = lrate*annealstep;          % anneal learning rate
            olddelta   = delta;                % accumulate angledelta until
            oldchange  = change;               %  annealdeg is reached
        elseif step == 1                     % on first step only
            olddelta   = delta;                % initialize
            oldchange  = change;
        end
        %
        %%%%%%%%%%%%%%%%%%%% Apply stopping rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if step >2 && change < nochange,      % apply stopping rule
            laststep=step;
            step=maxsteps;                  % stop when weights stabilize
        elseif change > DEFAULT_BLOWUP,      % if weights blow up,
            lrate=lrate*DEFAULT_BLOWUP_FAC;    % keep trying
        end;                                 % with a smaller learning rate
    end; % end if weights in bounds

end; % end training %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~laststep
    laststep = step;
end
lrates = lrates(1,1:laststep);           % truncate lrate history vector
%
%%%%%%%%%%%%%% Orient components towards max positive activation %%%%%%
%
if strcmp(posactflag,'on')
    [activations,winvout,weights] = posact(data,weights);
    % changes signs of activations and weights to make activations
    % net rms-positive
else
    activations = weights*data;
end
%
%%%%%%%%%%%%%% If pcaflag, compose PCA and ICA matrices %%%%%%%%%%%%%%%
%
if strcmp(pcaflag,'on')
    fprintf('Composing the eigenvector, weights, and sphere matrices\n');
    fprintf('  into a single rectangular weights matrix; sphere=eye(%d)\n'...
        ,chans);
    weights= weights*sphere*eigenvectors(:,1:ncomps)';
    sphere = eye(urchans);
end
%
%%%%%% Sort components in descending order of max projected variance %%%%
%
if verbose,
    fprintf(...
        'Sorting components in descending order of mean projected variance ...\n');
end
%
%%%%%%%%%%%%%%%%%%%% Find mean variances %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
meanvar  = zeros(ncomps,1);      % size of the projections
if ncomps == urchans % if weights are square . . .
    winv = inv(weights*sphere);
else
    fprintf('Using pseudo-inverse of weight matrix to rank order component projections.\n');
    winv = pinv(weights*sphere);
end
for s=1:ncomps
    if verbose,
        fprintf('%d ',s);         % construct single-component data matrix
    end
    % project to scalp, then add row means
    compproj = winv(:,s)*activations(s,:);
    meanvar(s) = mean(sum(compproj.*compproj)/(size(compproj,1)-1));
    % compute mean variance
end                                   % at all scalp channels
if verbose,
    fprintf('\n');
end
%
%%%%%%%%%%%%%% Sort components by mean variance %%%%%%%%%%%%%%%%%%%%%%%%
%
[sortvar, windex] = sort(meanvar);
windex = windex(ncomps:-1:1); % order large to small
meanvar = meanvar(windex);
%
%%%%%%%%%%%%%%%%%%%%% Filter data using final weights %%%%%%%%%%%%%%%%%%
%
if nargout>6, % if activations are to be returned
    if verbose,
        fprintf('Permuting the activation wave forms ...\n');
    end
    activations = activations(windex,:);
else
    clear activations
end
weights = weights(windex,:);% reorder the weight matrix
bias  = bias(windex);       % reorder them
signs = diag(signs);        % vectorize the signs matrix
signs = signs(windex);      % reorder them

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

return

%
%%%%%%%%%%%%%%%%%% return nonlinearly-transformed data  %%%%%%%%%%%%%%%%
%
if nargout > 7
    u=weights*data + bias*ones(1,frames);
    y = zeros(size(u));
    for c=1:chans
        for f=1:frames
            y(c,f) = 1/(1+exp(-u(c,f)));
        end
    end
end
