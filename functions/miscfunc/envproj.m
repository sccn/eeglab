% envproj() - plot envelopes of projections of selected ICA component 
%             projections against envelope of the original data
%
% Usage:    >> [envdata] = envproj(data,weights,compnums);
%           >> [envdata] = envproj(data,weights,compnums, ...
%                                  title,limits,chanlist,compnames,colors);
% Inputs:
%   data      = runica() input data (chans,frames) <- best one epoch only!
%   weights   = unmixing weight matrix (runica() weights*sphere)
%   compnums  = list of component numbers to project and plot
%
% Optional inputs:
%   title     = 'fairly short plot title' (in quotes) (0 -> none)
%   limits    = [xmin xmax ymin ymax] (x's in msec) 
%                       (0 or both y's 0 -> [min max])
%   chanlist  = list of data channels to use for max/min (0 -> all)
%   compnames = file of component labels (exactly 5 chars per line, 
%               trailing '.' = space) (default|0 -> compnums)
%   colors    = file of color codes, 3 chars per line  ('.' = space)
%                       (0 -> default color order (black/white first))
%
% Output: 
%   envdata   = [2,(frames*length(compnums)+1)] envelopes of data, comps.
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 1/1997 
%
% See also: envtopo()

% Copyright (C) 01-23-97 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 03-19-97 use datamean instead of frames/baseframes, diag(cov()) for var() -sm
% 04-03-97 changed name to envproj() -sm
% 05-20-97 used sum-of-squares for var instead of diag() to allow long data 
%          read MAXPLOTDATACHANS from icadefs.m -sm
% 06-07-97 changed order of args to conform to runica -sm
% 06-18-97 read MAXENVPLOTCHANS instead of MAXPLOTDATACHANS -sm
% 07-23-97 dropped datamean from args; mean is distributed among components -sm
% 10-11-97 range of max/min for default ylimits over chanlist only! -sm
% 11-05-97 disallowed white lines unless axis color is white -sm & ch
% 11-13-97 rm'd errorcode variable to improve limits handling -sm
% 11-18-97 added legend to plot; varied line types -sm
% 11-30-97 added f=figure below to keep previous plot from being altered -sm
% 12-01-97 added thicker lines for data envelope -sm
% 12-08-97 debugged thicker lines for data envelope -sm
% 12-10-97 changed f=figure below to f=gcf to avoid matlab bug (gcf confusion) -sm
% 12-10-97 implemented compnames arg -sm
% 12-13-97 took out compnames=compnames'; bug -sm
% 02-06-98 added 'show-largest' feature to legend when numcomps > 7 -sm
% 02-09-98 fixed numcomps bug -sm
% 02-12-98 test for outofbounds compnums early -sm
% 04-28-98 made to work for rectangular weights -sm
% 06-05-98 made NEG_UP optional -sm
% 02-22-99 added FILL mode default when numcomps==1 -sm
% 03-14-99 made envelope finding code more efficient -sm
% 12-19-00 updated icaproj() args -sm
% 01-25-02 reformated help & license, added links -ad 

function [envdata] = envproj(data,weights,compnums,titl,limits,chanlist,compnames,colors)

FILL = 1; % 1;  % 1 = fill in component envelope unless ncomps>1
DATA_LINEWIDTH = 1.5; % 2
COMP_LINEWIDTH = 1.5; % 1
MAXENVPLOTCHANS = 256;

if nargin < 3,
    help envproj
    fprintf('envproj(): must have at least three arguments.\n');
    return
end

icadefs        % load ENVCOLORS & MAXENVPLOTCHANS default variables
MAX_LEG = 7;   % maximum number of component line types to label in legend
NEG_UP = 0;    % 1 = plot negative up

%
%%%%%%%%%%%% Substitute for omitted arguments %%%%%%%%%%%%%%%%%%%%%%%%
%
if nargin < 8,
    colors = 'envproj.col';
elseif colors(1)==0,
    colors = 'envproj.col';
end

[chans,frames] = size(data);

if nargin < 7
   compnames = 0;
end

if nargin < 6
    chanlist = 0;
end
if nargin < 5,
    limits = 0;
end
if nargin < 4,
    titl = 0;
end
%
%%%%%%%%%%%%%%%%%%%%%%%% Test data size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
[wr,wc]       = size(weights);
%
%%%%%%%%%%%%%%%% Substitute for zero arguments %%%%%%%%%%%%%%%%%%%%%%%%%
%
if chanlist == 0,
    chanlist = [1:chans];
end
if compnums == 0,
    compnums = [1:wr];
end
if size(compnums,1)>1,        % handle column of compnums !
    compnums = compnums';
end
numcomps = length(compnums);
if numcomps > MAXENVPLOTCHANS,
    fprintf(...
'envproj(): cannot plot more than %d channels of data at once.\n',...
                  MAXENVPLOTCHANS);
    return
end

if max(compnums)>wr
 fprintf('\nenvproj(): Component index %d out of bounds (1:%d).\n',...
                                             max(compnums),wr);
 return
end

if min(compnums)<1,
 fprintf('\nenvproj(): Component index %d out of bounds (1:%d).\n',...
                                             min(compnums),wr);
 return
end

if FILL>0 && numcomps == 1
   FILL = 1;
else
   FILL = 0;
end
%
%%%%%%%%%%%%%%%%%%%% Read and adjust limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  if limits==0,
    xmin=0;xmax=1;
    ymin=min(min(data(chanlist,:)));
    ymax=max(max(data(chanlist,:)));
  else
    if length(limits)~=4,
        fprintf( ...
 'envproj(): limits should be 0 or an array [xmin xmax ymin ymax].\n');
        return
      end
        if limits(1,1) == 0 && limits(1,2) ==0,
            xmin=0;
            xmax=0;
        else
            xmin = limits(1,1);
            xmax = limits(1,2);
        end
         if limits(1,3) == 0 && limits(1,4) ==0,
            ymin=0;
            ymax=0;
        else
            ymin = limits(1,3);
            ymax = limits(1,4);
        end
  end

  if xmax == 0 && xmin == 0,
    x = (0:1:frames-1);
    xmin = 0;
    xmax = frames-1;
  else
    dx = (xmax-xmin)/(frames-1);
    x=xmin*ones(1,frames)+dx*(0:frames-1); % construct x-values
  end

  if ymax == 0 && ymin == 0,
    ymax=max(max(data(chanlist,:)));
    ymin=min(min(data(chanlist,:)));
  end

%
%%%%%%%%%%%%%%%%%%%% Read the color names %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    if ~ischar(colors)
        fprintf('envproj(): color file name must be a string.\n');
        return
    end
    cid = fopen(colors,'r');
    if cid <3,
        fprintf('envproj(): cannot open file %s.\n',colors);
        return
    else
        colors = fscanf(cid,'%s',[3 MAXENVPLOTCHANS]);
        colors = colors';
        [r c] = size(colors);
        for i=1:r
            for j=1:c
                if colors(i,j)=='.',
                    colors(i,j)=' ';
                end
            end
        end
    end
[rr cc] = size(colors);

%
%%%%%%%%%%%% Compute projected data for single components %%%%%%%%%%%%%
%
projdata = [data(chanlist,:)];
envdata = [min(projdata); max(projdata)]; % begin with envelope of data
fprintf('envproj(): Projecting component(s) ');
maxvars = zeros(1,length(compnums));
n=1;
for c=compnums,            % for each component 
  fprintf('%d ',c);
  % [icaproj] = icaproj(data,weights,compindex); % new arg order 12/00
  proj = icaproj(data,weights,c); % let offsets be 
                                         % distributed among comps.
  [val,i] = max(sum(proj.*proj)); % find max variance
  maxvars(n) = val;
  proj = proj(chanlist,:);
  projdata = [projdata proj];
  envtmp = [ min(proj) ; max(proj)];
  envdata = [envdata envtmp];   
                 % append envelope of projected data sets onto envdata
  % Note: size(envdata) = [length(chanlist)  frames*(numcomps+1)]
  n = n+1;
end
fprintf('\n');
%
%%%%%%%%%%%%%%%%%%%%%%%% Make the plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~ischar(titl)
  if titl==0,
    titl = ' ';
  else
    fprintf('envproj(): titl unrecognized.\n');
    return
  end
end

set(gcf,'Color',BACKCOLOR); % set the background color
set(gca,'Color','none'); % set the axis color to the background color
if NEG_UP
  set(gca,'Ydir','reverse');
end
n = 1;
H = zeros(1,numcomps+1);

for c=1:numcomps+1;
  %
  % Plot envelopes
  %
  if n <= 6
    if n == 1              % thick lines, store H for legend()
      H(n)=plot(x,envdata(1,(c-1)*frames+1:c*frames)',...
                                     colors(c),'Linewidth',DATA_LINEWIDTH);
      hold on;
           plot(x,envdata(2,(c-1)*frames+1:c*frames)',...
                                     colors(c),'Linewidth',DATA_LINEWIDTH);
    else
       if FILL
            fillx = [x x(frames:-1:1)];
            filly = [envdata(1,(c-1)*frames+1:c*frames) ...
                            envdata(2,c*frames:-1:(c-1)*frames+1)];
            i = find(isnan(filly)==1);
            if ~isempty(i)
                fprintf('Replacing %d NaNs in envelope with 0s.\n',length(i)/2);
                filly(i)=0;
            end
            H(n) = fill(fillx,filly,colors(c));
       else
            H(n)= plot(x,envdata(1,(c-1)*frames+1:c*frames)',...
                                     colors(c),'Linewidth',COMP_LINEWIDTH);
             hold on;
            plot(x,envdata(2,(c-1)*frames+1:c*frames)',...
                                     colors(c),'Linewidth',COMP_LINEWIDTH);
       end
    end
  else
    H(n)= plot(x,envdata(1,(c-1)*frames+1:c*frames)',...
                                     colors(c),'LineStyle',':','LineWidth',DATA_LINEWIDTH);
    hold on;
          plot(x,envdata(2,(c-1)*frames+1:c*frames)',...
                                     colors(c),'LineStyle',':','LineWidth',DATA_LINEWIDTH);
  end
  n = n+1;
end
set(gca,'Color','none'); % set the axis color to the background color

l=xlabel('Time (ms)');      % xaxis label
set(l,'FontSize',14);

l=ylabel('Potential (uV)'); % yaxis label
set(l,'FontSize',14);

if xmax > 0  && xmin < 0
  plot([0 0],[ymin ymax],'k','linewidth',1.5); % plot vertical line at time 0
end

if ymax ~= 0 |ymin~= 0,
    axis([min(x),max(x),ymin,ymax]);
else
    axis([min(x),max(x),min(envdata(1,:)'),max(envdata(2,:)')]);
end
legtext = [ 'Data ' ];
if compnames(1) ~= 0
    pid = fopen(compnames,'r');
    if pid < 3
       fprintf('envproj(): file %s not found.\n',compnames);
       return
    end
    compnames = fscanf(pid,'%s',[5 numcomps]);
    compnames = [ ['Data.']' compnames];
    compnames = compnames';
    if size(compnames,1) ~= numcomps+1
      fprintf('envproj(): no. of compnames must equal no. of compnums.\n');
      return
    end
    for c=1:5 % remove padding with .'s
      for r = 1:size(compnames,1)
         if compnames(r,c) == '.'
           compnames(r,c) = ' ';
         end
      end
    end
    legtext = compnames;
else
  legtext = ['Data'];
  for c = compnums
      legtext = strvcat(legtext,int2str(c));
  end
end

if numcomps<MAX_LEG+1
  if FILL <= 0  % no legend in FILL mode
       l=legend(H,legtext,0);      % plot legend
  end
else
  fprintf('Finding largest components for legend...\n');
  [maxvars,i] = sort(maxvars);
  maxvars = maxvars(length(maxvars):-1:1); % sort large-to-small
  i= i(length(i):-1:1);
  H = [H(1) H(i+1)];
  H = H(1:MAX_LEG+1);
  legtext = [legtext(1,:);legtext(i+1,:)];
  legtext = legtext(1:MAX_LEG+1,:);
  legend(H,legtext,0); % MAX_LEG largest-peaked comps
end

set(gca,'FontSize',12);
fprintf('If desired, use mouse to move legend or >> delete(legend).\n');

%     
%%%%%%%%%%%%% Compute percentage of variance accounted for %%%%%%%%%%%%
%
sumdata = zeros(length(chanlist),frames);
for e=2:numcomps+1,              % sum the component activities
    sumdata = sumdata + projdata(:,(e-1)*frames+1:e*frames);  
end

sigssqr = sum(data(chanlist,:).*data(chanlist,:))/(length(chanlist)-1);
dif    = data(chanlist,:)-sumdata;
difssqr = sum(dif.*dif)/(length(chanlist)-1);
pvaf   = round(100.0*(1.0-difssqr/sigssqr)); % % variance accounted for
rtitl = ['(' int2str(pvaf) '%)'];
%titl = [titl '   ' rtitl];

t=title(titl);              % plot title
set(t,'FontSize',14);

