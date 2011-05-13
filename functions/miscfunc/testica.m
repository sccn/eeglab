% testica() - Test the runica() function's ability to separate synthetic sources. 
%             Use the input variables to estimate the (best) decomposition accuracy
%             for a given data set size.
% Usage:
%        >> testica(channels,frames);  % No return variable -> plot results
%        >> [testresult] = testica(channels,frames,sources,exppow,shape);
%                          % Return variable -> return results with no plots
% Inputs:
%   channels = number of simulated data channels {no default}
%   frames   = number of simulated time points {no default}
%   sources  = number of simulated quasi-independent sources {default: =channels}
%   exppow   = exponential power for scaling size of the sources (0->all equal)
%               {default: -0.05 -> Ex: 14 sources scaled between 1.0 and 0.24}
%   shape    = varies monotonically with kurtosis of the simulated sources 
%               {default: 1.2 -> source kurtosis near 1 (super-Gaussian>0)}
%
% Authors: Scott Makeig & Te-Won Lee, SCCN/INC/UCSD, La Jolla, 2-27-1997 
%
% See also: runica()

% Copyright (C) 2-27-97 Scott Makeig & Te-Won Lee, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 2-28-97 added source, shape and exppow parameters, kurtosis -sm
% 4-03-97 shortened name to testica() -sm
% 4-14-97 changed call to runica() to use new variable order -sm
% 4-16-97 prints max and min abs corr -sm
% 7-10-97 changed to newrunica(), added surf() plot -sm
% 7-30-97 altered runica() call to fit version 3.0 -sm
% 3-02-00 replaced idit() call with call to Benjamin Blankertz' eyeLike() -sm
% 3-08-00 added kurt and exppow plots, changed defaults, added plot labels -sm
% 01-25-02 reformated help & license, added links -ad 

function [testresult] = testica(channels,frames,sources,exppow,shape)

icadefs; % read BACKCOLOR

% Default runica() parameter values:

block = 0;   % default block size 
lrate = 0;   % default starting lrate 
adeg  = 0;   % default annealing threshold 
maxsteps   = 0; % default 
sphereflag = 'on'; % default yes, perform sphering
stop  = 0.000001;  % default stopping wchange

% Defaults:

default_chans = 31;
default_frames = 10000;
% default sources = channels
default_exppow = -0.05;
default_shape = 1.2;
plotflag = 1;
try, plotflag = ismatlab; catch, end;

if nargin<2
  help testica
  return
end
if nargin<5
  shape = default_shape;
end
if nargin<4
   exppow = default_exppow;
end
if nargin<3,
	sources = 0;
end
if nargin<2
 frames = 0;
end
if nargin < 1
 channels   = 0;
end
if frames == 0,
  frames = default_frames;
end
if channels == 0,
  channels = default_chans;
end
if sources == 0,
	sources = channels;
end

if sources < channels,
  fprintf('testica() - sources must be >= channels.\n');
  exit 1
end

% Generate artificial super-Gaussian sources:

fprintf('\n  Testing runica() using %d simulated sources.\n\n',sources);
fprintf('Computing %d simulated source activations of length %d ...\n', ...
                                      sources,frames);
fprintf('Simulated source strengths: %4.3f to %4.3f.\n', ...
                                      1.0, exp(exppow*(channels-1)));
exppowers = zeros(1,channels);
exppowers(1) = 1.0;
for s=1:sources
  exppowers(s) = exp(exppow*(s-1));
end

% Synthesize random source activations

super=randn(sources,frames).*(exppowers'*ones(1,frames));  
super=sign(super).*abs(super.^shape); % make super-Gaussian if shape > 1
% fprintf('Size of super = %d,%d\n',size(super,1),size(super,2));

if frames > 40 && plotflag
  figure
  pos = get(gcf,'position');
  off = [40 -40 0 0]; % succeeding figure screen position offsets

  hist(super(1,:),round(frames/20));
  tt=title('Amplitude distribution of source 1');
  set(tt,'fontsize',14);
  xlm=get(gca,'xlim');
  ylm=get(gca,'ylim');
  kurttext = ['Kurtosis = ' num2str(kurt(super(1,:)),3)];
  tp=[xlm;ylm]*[0.25;0.75];
  kt=text(tp(1),tp(2),kurttext);
  set(kt,'fontsize',13);
  set(kt,'horizontalalignment','center');
else
  fprintf('Not plotting source amplitude histogram: data length too small.\n')
end 

if nargout == 0 && plotflag
  input('Hit enter to view source strengths: ');
  fprintf('\n')
  if frames <= 40
    figure
    pos = get(gcf,'position');
  else
    figure('position',pos+off);
  end
  plot(1:sources,exppowers);
  hold on;plot(1:sources,exppowers,'r^');
  set(gca,'xlim',[0 sources+1]);
  set(gca,'ylim',[0 1]);
  xt=title(['Relative source amplitudes (exppow = ' num2str(exppow,3) ')']);
  set(xt,'fontsize',14);
  axl=xlabel('Source Number');
  ayl=ylabel('Relative Amplitude');
  set(axl,'fontsize',14);
  set(ayl,'fontsize',14);
end

k = kurt(super'); % find kurtosis of rows of super
maxkurt = max(k); minkurt=min(k);
fprintf('Simulated source kurtosis: %4.3f to %4.3f.\n',minkurt,maxkurt);

tmp = corrcoef(super');
i = find(tmp<1);
minoff = min(abs(tmp(i)));
maxoff = max(abs(tmp(i)));
fprintf('Absolute correlations between sources range from %5.4f to %5.4f\n', ...
                                                           minoff,maxoff);

fprintf('Mixing the simulated sources into %d channels ...\n',channels);
forward = randn(channels,sources);   % random forward mixing matrix
data = forward*super; % these are the simulated observed data
if nargout == 0
    input('Hit enter to start ICA decomposition: ')
    fprintf('\n')
end

fprintf('Decomposing the resulting simulated data using runica() ...\n');
 [weights,sphere,compvars,bias,signs,lrates,activations] = runica(data, ...
        'block',block, ...
           'lrate',lrate, ...
              'nochange',stop, ...
                 'annealdeg',adeg, ...
                    'maxsteps',maxsteps, ...
                       'sphering',sphereflag, 'weights',eye(channels)); 
 
fprintf('ICA decomposition complete.\n');

% Alternatively, activations = icaact(data,weights,sphere,datamean);

fprintf('\nScaling and row-permuting the resulting performance matrix ...\n')
testid = weights*sphere*forward(:,1:channels); 
                              % if separation were complete, this
                              % would be a scaled and row-permuted 
                              % identity matrix
testresult = eyelike(testid); % permute output matrix rows to rememble eye()
                              % using Benjamin Blankertz eyeLike() - 3/2/00
% testresult = idit(testid);  % permute output matrix rows to rememble eye()
                              % scale to make max column elements all = 1

tmp = corrcoef(activations');
i = find(tmp<1);
maxoff = max(abs(tmp(i)));
minoff = min(abs(tmp(i)));
fprintf('Absolute activation correlations between %5.4f and %5.4f\n',minoff,maxoff);

i = find(testresult<1);       % find maximum abs off-diagonal value
maxoff = max(abs(testresult(i)));
meanoff = mean(abs(testresult(i)));
if sources > channels,
 fprintf('The returned matrix measures the separation');
 fprintf('of the largest %d simulated sources,\n',channels);
end

fprintf('Perfect separation would return an identity matrix.\n');
fprintf('Max absolute off-diagonal value in the returned matrix: %f\n',maxoff);
fprintf('Mean absolute off-diagonal value in the returned matrix: %f\n',meanoff);

 [corr,indx,indy,corrs] = matcorr(activations,super);

fprintf('Absolute corrs between best-matching source and activation\n');
fprintf('         component pairs range from %5.4f to %5.4f\n', ...
                                 abs(corr(1)),abs(corr(length(corr))));

if nargout == 0
   fprintf('\nView the results:\n');
   fprintf('Use mouse to rotate the image.\n');
end
if ~plotflag, return; end;
figure('Position',pos+2*off);
set(gcf,'Color',BACKCOLOR);
surf(testresult);  % plot the resulting ~identity matrix
st=title('Results: Test of ICA Separation');
set(st,'fontsize',14)
sxl=xlabel('Source Out');
set(sxl,'fontsize',14);
syl=ylabel('Source In');
set(syl,'fontsize',14);
szl=zlabel('Relative Recovery');
set(szl,'fontsize',14);
view(-52,50);
axis('auto');
if max(max(abs(testresult)))>1
  fprintf('NOTE: Some sources not recovered well.\n');
  fprintf('Restricting plot z-limits to: [-1,1]\n');
  set(gca,'zlim',[-1 1]);
end
rotate3d

if nargout == 0
   testresult = [];
end
