% compdsp() - Display standard info figures for a data decomposition
%             Creates four figure windows showing: Component amplitudes,
%             scalp maps, activations and activation spectra.
% Usage:
%      >> compdsp(data,weights,locfile,[srate],[title],[compnums],[amps],[act]);
%
% Inputs:
%   data     = data matrix used to train the decomposition
%   weights  = the unmixing matrix (e.g., weights*sphere from runica())
%
% Optional:
%   locfile  = 2-d electrode locations file (as in >> topoplot example)
%              {default|[]: default locations file given in icadefs.m
%   srate    = sampling rate in Hz {default|[]: as in icadefs.m}
%   title    = optional figure title text 
%              {default|'': none}
%   compnums = optional vector of component numbers to display 
%              {default|[] -> all}
%   amps     = all component amplitudes (from runica()) 
%              {default|[]->recompute}
%   act      = activations matrix (from runica())
%              {default|[]->recompute}
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 2000 

% Copyright (C) 12/16/00 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 02-01-01  replaced std() with rms() -sm
% 02-10-01  made srate optional -sm
% 01-25-02 reformated help & license -ad 

function compdsp(data,unmix,locfile,srate,titl,compnums,amps,act)

minHz = 2; % frequency display limits
maxHz = 40;
ACTS_PER_EEGPLOT = 32;

%
%%%%%%%%%%%%%%%%%%%%%% Read and test arguments %%%%%%%%%%%%%%%%%%%%%%%%
%
icadefs % read BACKCOLOR, DEFAULT_SRATE

if nargin<2
  help compdsp
  return
end

chans = size(data,1);
frames = size(data,2);
ncomps = size(unmix,1);

if ncomps < 1
   error('Unmixing matrix must have at least one component');
end
if chans < 2
   error('Data must have at least two channels');
end
if frames < 2
   error('Data must have at least two frames');
end
if size(unmix,2) ~= chans
   error('Sizes of unmix and data do not match');
end
if nargin<3
 locfile=[];
end
if isempty(locfile)
   locfile = DEFAULT_ELOC; % from icsdefs.m
end
if ~exist(locfile)
   error('Cannot find electrode locations file');
end

if nargin<4
  srate = 0; % from icadefs
end
if isempty(srate) | srate==0
  srate = DEFAULT_SRATE; % from icadefs
end
if nargin<5
  titl = '';    % default - no title text
end
if isempty(titl)
  titl = '';
end
if nargin<6
  compnums = 0; % default - all components
end
if isempty(compnums)
  compnums = 0;
end
if nargin<7
  amps = NaN; % default - recompute amps
end
if isempty(amps)
  amps = NaN;
end
if ~isnan(amps) & length(amps) ~= ncomps
   error('Supplied amps does not match the number of unmixed components');
end

if compnums(1) == 0
  compnums = 1:ncomps;
end
if min(compnums) < 1
   error('Compnums must be positive');
end
if min(compnums) > ncomps
   error('Some compnum > number of components in unmix');
end

if nargin<8
   act = NaN; % default - recompute act
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Create plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if chans == ncomps 
  winv = inv(unmix);
elseif chans < ncomps 
   error('More components than channels?');
else
  winv = pinv(unmix);
end

if isnan(act)
  fprintf('Computing activations ...')
  act = unmix(compnums,:)*data;
  fprintf('\n');
elseif size(act,2) ~= frames
  error('Supplied activations do not match data length');
elseif size(act,1) ~=  ncomps & size(act,1) ~= length(compnums)
  error('Number of supplied activations matrix does not match data or weights');
elseif size(act,1) ==  ncomps 
  act = act(compnums,:); % cut down matrix to desired components
end

%
%%%%%%%%%%%%%%%%%%%%% I. Plot component amps %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
pos = [40,520,550,400]; figure('Position',pos);
if isnan(amps)
  fprintf('Computing component rms amplitudes ');
  amps = zeros(1,length(compnums));
  for j=1:length(compnums)
     amps(j) = rms(winv(:,compnums(j)))*rms(act(j,:)');
     fprintf('.')
  end
  fprintf('\n');
else
  amps = amps(compnums); % truncate passed amps to desired compnums
end
plot(compnums,amps,'k');
hold on
plot(compnums,amps,'r^');
xl=xlabel('Component numbers');
yl=ylabel('RMS Amplitudes');
tl=title([titl ' Amplitudes']);
ax = axis;
axis([min(compnums)-1 max(compnums)+1 0 ax(4)]);

%
%%%%%%%%%%%%%%%%%%%%%%%%% II. Plot component maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
pos = [40,40,550,400]; figure('Position',pos);

fprintf('Plotting component scalp maps ...') % compmaps() may make multiple figures
compmap(winv,locfile,compnums,[titl ' Scalp Maps'],0,compnums); 
fprintf('\n');

%
%%%%%%%%%%%%%%%%%%%%%% III. eegplot() activations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if frames/srate < 10
  dispsecs = ceil(frames/srate);
else
  dispsecs = 10; % defaults - display 10s data per screen
end
range = 0.8*max(max(act')-min(act'));

stact=1;
lact=ACTS_PER_EEGPLOT;
if lact>size(act,1)
   lact = size(act,1);
end

pos = [620,520,550,400]; figure('Position',pos);
while stact <= size(act,1)
  % eegplot(data,srate,spacing,eloc_file,windowlength,title,posn)
  eegplot(act(stact:lact,:),srate,range,compnums(stact:lact),...
               dispsecs,[titl ' Activations'],pos);
  pos = pos + [.02 .02 0 0];
  stact = stact+ACTS_PER_EEGPLOT;
  lact  = lact +ACTS_PER_EEGPLOT;
  if lact>size(act,1)
   lact = size(act,1);
  end
end

%
%%%%%%%%%%%%%%%%%%%%%%% IV. plotdata() spectra %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
pos = [620,40,550,400]; figure('Position',pos);

if frames > 2048
  windw = 512;
elseif frames > 1024
  windw = 256;
else
  windw = 128;
end
fprintf('Computing component spectra ')
for j = 1:length(compnums)
  % [Pxx,F] = PSD(X,NFFT,Fs,WINDOW,NOVERLAP)
  [spec,freqs] = psd(act(j,:),1024,srate,windw,ceil(windw*0.5));
  if ~exist('specs')
     specs = zeros(length(compnums),length(freqs));
  end
  specs(j,:) = spec';
  fprintf('.')
end
fprintf('\n');
specs = 10*log10(specs);

tmp = ceil(sqrt(length(compnums)));
tmp2 = ceil(length(compnums)/tmp);
for j=1:length(compnums)
  sbplot(tmp2,tmp,j)
  plot(freqs,specs(j,:))
  set(gca,'box','off')
  set(gca,'color',BACKCOLOR);
  ax=axis;
  axis([minHz,maxHz,ax(3),ax(4)]);
  tl = title(int2str(compnums(j)));
end
xl=xlabel('Frequency (Hz)');
yl=ylabel('Power (dB)');
set(gca,'YAxisLocation','right');
txl=textsc([titl ' Activation Spectra'],'title');
axcopy % pop-up axes on mouse click

% plottopo(specs,[tmp2 tmp],0,[2,70 min(min(specs)) max(max(specs))],...
%                 [titl ' Activation Power Spectra']);

function rmsval = rms(column)
   rmsval = sqrt(mean(column.*column)); % don't remove mean
