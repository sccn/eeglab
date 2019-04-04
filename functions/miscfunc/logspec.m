% logspec() - plot mean log power spectra of submitted data on loglog scale
%             using plotdata() or plottopo() formats 
%
% Usage:
%        >> [spectra,freqs] = logspec(data,frames,srate);
%        >> [spectra,freqs] = logspec(data,frames,srate,'title',...
%                                    [loHz-hiHz],'chan_locs',rm_mean);
% Inputs:
%    data   = input data (chans,frames*epochs)
%    frames = data samples per epoch   {default length(data)}
%    srate  = data sampling rate in Hz {default 256 Hz};
%    'title' = plot title {default: none}
%    [loHz-hiHz] = [loHz hiHz] plotting limits 
%       {default: [srate/fftlength srate/2]}
%    'chan_locs' = channel location file (ala topoplot()) 
%                   Else [rows cols] to plot data in a grid array
%    rm_mean = [0/1] 1 -> remove log mean spectrum from all
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 11-07-97 
%
% See also: plotdata(), plottopo()

% Copyright (C) 11-07-97 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% Changed plotdata() below to plottopo() 11/12/99 -sm
% Mentioned new grid array option 12/22/99 -sm
% 01-25-02 reformated help & license, added links -ad 

function [spectra,freqs] = logspec(data,frames,srate,titl,Hzlimits,chanlocs,rmmean)

if nargin < 1
  help logspec
  return
end
[rows,cols] = size(data);

icadefs  % read plotdata chan limit

if rows > MAXPLOTDATACHANS
   fprintf('logspec(): max plotdata() channels is %d.\n',MAXPLOTDATACHANS);
   return
end
if rows < 2
   fprintf('logspec(): min plotdata() channels is %d.\n',2);
   return
end
if nargin<7
  rmmean = 0; % default
end

if nargin < 5
 Hzlimits = 0;
end

if nargin < 4
  titl = ' ';
end

if nargin < 3
  srate = 256;
end

if nargin < 2,
  frames = cols;
end

epochs = fix(cols/frames);
if epochs*frames ~= cols
    fprintf('logspec() - frames does not divide data length.\n');
    return
end

fftlength = 2^floor(log(frames)/log(2));

spectra = zeros(rows,epochs*fftlength/2);
f2 = fftlength/2;
dB = 10/log(10);

if length(Hzlimits) < 2
 Hzlimits = [srate/fftlength srate/2];
end
if Hzlimits(2) <= Hzlimits(1)
   help logspec
   return
end

for e=1:epochs
  for r=1:rows
    [Pxx,freqs] = pwelch(data(r,(e-1)*frames+1:e*frames),fftlength,...
                          fftlength/4,fftlength,srate);
    spectra(r,(e-1)*f2+1:e*f2) = Pxx(2:f2+1)'; % omit DC bin
  end
end

clf
freqs = freqs(2:f2+1);
fsi = find(freqs >= Hzlimits(1) & freqs <= Hzlimits(2));
minf = freqs(fsi(1));
maxf = freqs(fsi(length(fsi)));
nfs = length(fsi);

showspec = zeros(rows,length(fsi)*epochs);
for e = 1:epochs
 showspec(:,(e-1)*nfs+1:e*nfs) = dB*log(spectra(:,(e-1)*f2+fsi));
end
% minspec = min(min(showspec));
% showspec = showspec-minspec; % make minimum 0 dB

showspec = blockave(showspec,nfs);

% meanspec = mean(showspec);
% showspec = showspec - ones(rows,1)*meanspec;

% >> plotdata(data,frames,limits,title,channames,colors,rtitle)
% diff = 0;
% MINUEND = 6;
% for r=1:rows
%  diff = diff - MINUEND;
%  showspec(r,:) = showspec(r<:)-diff;
% end
% semilogx(freqs(fsi),showspec');
% ax = axis;
% axis([minf maxf ax(3) ax(4)]);
% title(titl);

if nargin<6
% >> plotdata(data,frames,limits,title,channames,colors,rtitle,ydir)
     if rmmean
        showspec = showspec - ones(rows,1)*mean(showspec);
     end
     plotdata(showspec,nfs,[minf maxf 0 0],titl);
else
% >> plottopo(data,'chan_locs',frames,limits,title,channels,axsize,colors,ydir)
     if rmmean
        showspec = showspec - ones(rows,1)*mean(showspec);
     end
     plottopo(showspec,chanlocs,nfs,[minf maxf 0 0],titl);
end
ax = get(gcf,'children');
for a = ax
  set(a,'XScale','log')
end
