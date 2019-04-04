% chanproj() - make a detailed plot of data returned from plotproj() 
%              for given channel. Returns the data plotted.
% Usage:
%  >> [chandata] = chanproj(projdata,chan);
%  >> [chandata] = chanproj(projdata,chan,ncomps,framelist,limits,title,colors);
%
% Inputs:
%   projdata    = data returned from plotproj() 
%   chan        = single channel to plot
%   ncomps      = number of component projections in projdata
%
% Optional:
%   framelist   = data frames to plot per epoch Ex: [1:128] (0|def -> all)
%   limits      = [xmin xmax ymin ymax]  (x's in msec) 
%                                        (0, or y's 0 -> data limits)
%   title       = fairly short single-quoted 'plot title'  (0|def -> chan)
%   colors      = file of color codes, 3 chars per line  (NB: '.' = space) 
%                                        (0|def -> white is original data)
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 1996 

% Copyright (C) 1996 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 11-30-96 Scott Makeig  CNL / Salk Institute, La Jolla as plotprojchan.m
% 03-19-97 changed var() to diag(cov()) -sm
% 03-26-97 fix (== -> =) -sm
% 04-03-97 allow framelist to be a col vector, allow 32 traces, fix pvaf, made
%           ncomps mandatory  -sm
% 04-04-97 shortened name to chanproj() -sm
% 05-20-97 added read of icadefs.m -sm
% 11-05-97 disallowed white traces unless default axis color is white -sm & ch
% 11-13-97 rm'ed errcode variable -sm
% 12-08-97 added LineWidth 2 to data trace, changed whole plot color to BACKCOLOR -sm
% 01-25-02 reformated help & license -ad 

function [chandata] = chanproj(projdata,chan,ncomps,framelist,limits,titl,colorfile)

icadefs; % read default MAXPLOTDATACHANS

if nargin < 7,
    colorfile =0;
end
if nargin < 6,
    titl = 0;
end
if nargin < 5,
    limits = 0;
end
if nargin < 4,
    framelist = 0;
end

if nargin < 3,
    fprintf('chanproj(): requires at least three arguments.\n\n');
    help chanproj
    return
end

[chans,framestot] = size(projdata);
frames = fix(framestot/(ncomps+1));

if ncomps < 1 || frames*(ncomps+1) ~= framestot,
    fprintf(...
  'chanproj(): data length (%d) not a multiple of ncomps (%d).\n',...
                                   framestot,ncomps);
    return
end
 
if chan> chans,
    fprintf(...
 'chanproj(): specified channel (%d) cannot be > than nchans (%d).\n',...
                            chan,chans);
    return
else
    chandata = projdata(chan,:);
    epochs = fix(length(chandata)/frames);
end
if epochs > MAXPLOTDATACHANS
  fprintf(...
 'chanproj(): maximum number of traces to plot is %d\n',...
              MAXPLOTDATACHANS);
   return
end
if framestot~=epochs*frames,
   fprintf('chanproj(): projdata is wrong length.\n');  % exit
   return
end

%
%%%%%%%%%%%%%%%% Find or read plot limits %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if limits==0,
  xmin=0;xmax=0;ymin=0;ymax=0;
else
  if length(limits)~=4,
      fprintf( ...
'chanproj():^G limits should be 0 or an array [xmin xmax ymin ymax].\n');
      return
  end
  xmin = limits(1);
  xmax = limits(2);
  ymin = limits(3);
  ymax = limits(4);
end

if xmin == 0 && xmax == 0,
  x = [0:frames-1];
  xmin = 0;
  xmax = frames-1;
else
  dx = (xmax-xmin)/(frames-1);
  x=xmin*ones(1,frames)+dx*(0:frames-1);          % construct x-values
end

if ymax == 0 && ymin == 0,
  ymax=max(max(projdata(chan,:)));
  ymin=min(min(projdata(chan,:)));
end
%
%%%%% Reshape the projdata for plotting and returning to user %%%%%%%%%
%
chandata=reshape(chandata,frames,epochs);
chandata=chandata';

if framelist ==0,
  framelist = [1:frames];  % default
end
if size(framelist,2)==1,
     if size(framelist,1)>1,
         framelist = framelist';
     else
         fprintf(...
 'chanproj(): framelist should be a vector of frames to plot per epoch.\n');
         return
     end
end
  
if framelist(1)<1 || framelist(length(framelist))>frames,
    fprintf('chanproj(): framelist values must be between 1-%d.\n',frames);
    return
else
    chandata=chandata(:,framelist);
end
%
%%%%%%%%%%%%%%%%%%% Read the color names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if colorfile ~=0,
    if ~ischar(colorfile)
       fprintf('chanproj(); color file name must be a string.\n');
       return
    end
    cid = fopen(colorfile,'r');
    if cid <3,
        fprintf('chanproj(): cannot open file %s.\n',colorfile);
        return
    end
    colors = fscanf(cid,'%s',[3 MAXPLOTDATACHANS]);
    colors = colors';
    [r c] = size(colors);
    for i=1:r
        for j=1:c
            if colors(i,j)=='.',
                colors(i,j)=' ';
            end
        end
    end
else    % default color order - no yellow
    colors =['w  ';'r  ';'b  ';'g  ';'c  ';'m  ';'r  ';'b  ';'g  ';'c  ';'m  ';'r  ';'b  ';'g  ';'c  ';'m  ';'r  ';'b  ';'g  ';'c  ';'m  ';'r  ';'b  ';'g  ';'c  ';'m  ';'r  ';'b  ';'g  ';'c  ';'m  ';'r  '];
end
%
% Make vector of x-values
%
x=[xmin:(xmax-xmin)/(frames-1):xmax+0.00001];
xmin=x(framelist(1));
xmax=x(framelist(length(framelist)));
x = x(framelist(1):framelist(length(framelist)));
%
%%%%%%%%% Compute percentage of variance accounted for %%%%%%%%%%%%%%
%
if epochs>1,
    sumdata = zeros(1,length(framelist));
    for e=2:epochs
        sumdata= sumdata + chandata(e,:);  % sum the component projections
    end
    sigvar = diag(cov(chandata(1,:)'));
    difvar = diag(cov(((chandata(1,:)-sumdata)')));
                                    % percent variance accounted for
    pvaf   = round(100.0*(1.0-difvar/sigvar)); 
end
%
%%%%%%%%%%%%%%%%%%%%%% Plot traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
fprintf('chanproj(): Drawing trace ');
for e=1:epochs,
    fprintf ('%d ',e);
    set(gcf,'Color',BACKCOLOR); % set the background color to grey
    set(gca,'Color','none');    % set the axis color = figure color
    if e==1
      plot(x,chandata(e,:),colors(e),'LineWidth',2);  % plot it!
    else
      plot(x,chandata(e,:),colors(e),'LineWidth',1);  % plot it!
    end
    hold on;
end
fprintf('\n');
%
%%%%%%%%%%%%%%%%%%%% Fix axis limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
axis([xmin xmax ymin-0.1*(ymax-ymin) ymax+0.1*(ymax-ymin)]);
%
%%%%%%%%%%%%%%%%%%% Add title and labels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if titl==0,
    titl = ['channel ' int2str(chan)];
end
titl = [ titl '   (p.v.a.f. ' int2str(pvaf) '%)' ];
title(titl);
xlabel('Time (msec)');
ylabel('Potential (uV)');
