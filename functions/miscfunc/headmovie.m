% ########## This function is deprecated. Use eegmovie instead. ########## 
%
% headmovie() - Record a Matlab movie of scalp data.
%               Use seemovie() to display the movie.
%
% Usage:  >> [Movie,Colormap] = headmovie(data,elec_loc,spline_file);
%         >> [Movie,Colormap,minc,maxc] = headmovie(data,elec_loc,spline_file,...
%                                srate,title,camerapath,movieframes,minmax,startsec);
% Inputs:
%   data        = (chans,frames) EEG data set to plot
%   elec_loc    = electrode locations file for eegplot {default 'chan.loc'}
%   spline_file = headplot() produced spline 'filename' {default  'chan.spline'}
%   srate       = sampling rate in Hz {default|0 -> 256 Hz}
%   title       = 'plot title' {default|0 -> none}
%   camerapath  = [az_start az_step el_start el_step] {default [-127 0 30 0]}
%                 Setting all four non-0 creates a spiral camera path
%                 Partial entries allowed, e.g. [az_start]
%                 Subsequent rows [movieframe az_step 0 el_step] adapt step 
%                 sizes allowing starts/stops, panning back and forth, etc.
%   movieframes = vector of data frames to animate {default|0 -> all}
%   minmax      = Data [lower_bound, upper_bound] Lower->blue, upper->red
%                 {default|0 -> +/-abs max of data}
%   startsec    = starting time in seconds {default|0 -> 0.0}
%
% Note: BUG IN MATLAB 5.0-5.1 -- CANT SHOW MOVIES IN CORRECT COLORS
%
% Authors: Scott Makeig & Colin Humphries, SCCN/INC/UCSD, La Jolla, 2/1998 
%
% See also: seemovie(), eegmovie(), headplot()

% Copyright (C) 2.6.98 Scott Makeig & Colin Humphries, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 01-25-02 reformated help & license, added links -ad 

function [Movie, Colormap, minc, maxc] = headmovie(data,eloc_file,spline_file,srate,titl,camerapath,movieframes,minmax,startsec,varargin)

if nargin<1
	help headmovie
    return
end

clf
[chans,frames] = size(data);

icadefs;   % read DEFAULT_SRATE;

if nargin<9
   startsec = 0;
end
if nargin<8
   minmax = 0;
end
if size(minmax)==[1 1] && minmax ~= 0
   minmax = [-minmax minmax]; 
end
if minmax ==0,
	minc = min(min(data));
	maxc = max(max(data));
    absmax  = max([abs(minc), abs(maxc)]);
    fudge   = 0.05*(maxc-minc); % allow for slight extrapolation
    minc = -absmax-fudge;
    maxc =  absmax+fudge;
    minmax  = [minc maxc];
end

if nargin <7
	movieframes = 0;
end
if nargin<6
    camerapath = 0;
end
if movieframes == 0
	movieframes = 1:frames; % default
end
if nargin <5
	titl = '';
end
if titl == 0
	titl = '';
end
if nargin <4
	srate = 0;
end
if nargin <3
	spline_file = 0;
end
if nargin <2
	eloc_file = 0;
end

if movieframes(1) < 1 || movieframes(length(movieframes))>frames
	fprintf('headmovie(): specified movieframes not in data!\n');
	return
end

if srate ==0,
	srate = DEFAULT_SRATE;
end

mframes = length(movieframes);
fprintf('Making a movie of %d frames\n',mframes)
Movie    = moviein(mframes,gcf);

if size(camerapath) == [1,1]
   if camerapath==0
     camerapath = [-127 0 30 0];
     fprintf('Using default view [-127 0 30 0].');
   else
     camerapath = [camerapath 0 30 0];
   end
elseif size(camerapath,2) == 2
     camerapath(1,:) = [camerapath(1,1) camerapath(1,2) 30 0]
elseif size(camerapath,2) == 3
     camerapath(1,:) = [camerapath(1,1) camerapath(1,2) camerapath(1,3) 0]
elseif size(camerapath,2) > 4
     help headmovie
     return
end
if size(camerapath,1) > 1 && size(camerapath,2)~=4
     help headmovie
     return
end

azimuth   = camerapath(1,1); % initial camerapath variables
az_step   = camerapath(1,2);
elevation = camerapath(1,3);
el_step   = camerapath(1,4);

if ~isstruct(eloc_file) && eloc_file(1) == 0
   eloc_file    = 'chan.angles';
end
if spline_file(1) == 0
   spline_file  = 'chan.spline';
end

figure(gcf);   % bring figure to front
clf            % clear figure

%%%%%%%%%%%%%%%%%%%%% eegplot() of data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axeegplot = axes('Units','Normalized','Position',[.75 .05 .2 .9]);
if isstruct(eloc_file)
    eegplotold('noui',-data,srate,0,[],startsec,'r');
    %eegplotsold(-data,srate,[],' ',0,frames/srate,'r',startsec); %CJH
else
    eegplotold('noui',-data,srate,0,eloc_file,startsec,'r');
    %eegplotsold(-data,srate,eloc_file,' ',0,frames/srate,'r',startsec); %CJH
end;    
set(axeegplot,'XTick',[])                %%CJH
                                         % plot negative up
limits = get(axeegplot,'Ylim');          % list channel numbers only
set(axeegplot,'GridLineStyle',':')
set(axeegplot,'Xgrid','off')
set(axeegplot,'Ygrid','on')
axcolor = get(gcf,'Color');

%%%%%%%%%%%%%%%%%%%%%%% print title text %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axtitle    = axes('Units','Normalized','Position',[0 .9 .72 .1],'Color',axcolor);
   axes(axtitle)
   text(0.5,0.5,titl,'FontSize',16,'horizontalalignment','center')
   axis('off')

%%%%%%%%%%%%%%%%%%%%%%% display colorbar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axcolorbar = axes('Units','Normalized','Position',[.05 .05 .055 .18]); 
   axes(axcolorbar)
   if exist('colorbar_tp')==2
     h=colorbar_tp(axcolorbar); % Slightly altered Matlab colorbar() routine
   else
     h=colorbar(axcolorbar);  % Note: there is a minor problem with this call.
   end                 % Write authors for further information.
   hkids=get(h,'children');
   for k=1:1
      delete(hkids(k));
   end
   set(axcolorbar,'Ytick',[0 0.5 1.0],'Yticklabel',[' - ';' 0 ';' + '],'FontSize',16);
   axis('off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%% "Roll'em!" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axheadplot = axes('Units','Normalized','Position',[0 .1 .72 .8],'Color',axcolor);
TPAXCOLOR  = get(axheadplot,'Color');    
Colormap   = [jet(64);TPAXCOLOR];       

fprintf('Warning: do not resize plot window during movie creation ...\n   ');

newpan = length(movieframes)+1;
posrow = 2;
if size(camerapath,1) > 1
  newpan = camerapath(posrow,1); % pick up next frame to change camerapos step values
end
for i = movieframes                      % make the movie, frame by frame
   fprintf('%d ',i-movieframes(1)+1);

   axes(axeegplot)
   x1 = startsec+(i-1)/srate;
   l1 = line([x1 x1],limits,'color','b'); % draw vertical line at data frame
   timetext = num2str(x1,3);   % Note: Matlab4 doesn't take '%4.3f'
   set(axeegplot,'Xtick',x1,'XtickLabel',num2str(x1,'%4.3f'));

   axes(axheadplot)
   cla
   set(axheadplot,'Color',axcolor);

   headplot(data(:,i),spline_file,'maplimits',minmax,...
                  'view',[azimuth elevation],'verbose','off', varargin{:}); 

   if i== newpan    % adapt camerapath step sizes 
     az_step = camerapath(posrow,2);
     el_step = camerapath(posrow,4);
     posrow = posrow+1;
     if size(camerapath,1)>=posrow
       newpan = camerapath(posrow,1);
     else
       newpan = length(movieframes)+1;
     end
   end

   azimuth = azimuth+az_step;     % update camera position
   elevation = elevation+el_step;

   if elevation>=90 
      fprintf('headplot(): warning -- elevation out of range!\n');
      elevation = 89.99;
   end
   if elevation<=-90 
      fprintf('headplot(): warning -- elevation out of range!\n');
      elevation = -89.99;
   end

   set(axheadplot,'Units','pixels',...
    'CameraViewAngleMode','manual',... 
    'YTickMode','manual','ZTickMode','manual',...
    'PlotBoxAspectRatioMode','manual',...
    'DataAspectRatioMode','manual');    % keep camera distance constant

   txt = [int2str(i-movieframes(1)+1)];  
   text(0.45,-0.45,txt,'Color','k','FontSize',14);    % show frame number

   [Movie(:,i+1-movieframes(1))] = getframe(gcf); % can't show this - colors wrong!
                                                  % known Matlab 5.0-5.1 bug!
   % [M,Cmap] = getframe(gcf);
   % if i==1
   %   [m,n] = size(M);
   %   Movie = zeros(m,n*length(movieframes));
   %   whos Movie
   % end
   % Movie((i-1)*m+1:i*m,(i-1)*n+1:i*n) = M;
   drawnow
   delete(l1)
end
% colormap(Cmap);
% montage(Movie)
% Movie = immovie(Movie,[m n length(movieframes)]);
%    whos Movie
fprintf('\nDone!\n');

