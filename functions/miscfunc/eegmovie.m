% eegmovie() - Compile and view a Matlab movie. 
%              Uses scripts eegplotold() and topoplot().
%              Use seemovie() to display the movie.
% Usage:
%
% >> [Movie,Colormap] = eegmovie(data,srate,elec_locs,title,movieframes,minmax,startsec,...);
%
% Inputs:
%   data        = (chans,frames) EEG data set to plot
%   srate       = sampling rate in Hz {0 -> 256 Hz}
%   elec_locs   = ascii file of electrode locations {0 -> 'chan_file'}
%   title       = 'plot title' {0 -> none}
%   movieframes = vector of frames to animate {0 -> all}
%   minmax      = [blue_lower_bound, red_upper_bound] 
%                 {0 -> +/-abs max of data}
%  startsec     = starting time in seconds {0 -> 0.0}
%  additional options from topoplot are allowed
%  
%
% Author: Colin Humphries & Scott Makeig, CNL, Salk Institute, La Jolla, 3/97
%
% See also: seemovie(), eegplotold(), topoplot()

% Copyright (C)  6/4/97 Colin Humphries & Scott Makeig, CNL / Salk Institute / La Jolla CA
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

% 6/6/97 added movieframes arg -sm
% 6/12/97 removed old 'startframes' var., fixed vertical line frame selection -sm
% 6/27/97 debugged vertical line position -sm
% 10/4/97 clarified order of srate and eloc_locs -sm
% 3/18/97 changed eegplots -> eegplot('noui') -sm
% 10/10/99 added newlines to frame print at suggestion of Ian Lee, Singapore -sm
% 01/05/01 debugged plot details -sm
% 01/24/02 updated eegplot to eegplotold -ad
% 01-25-02 reformated help & license, added links -ad 

function [Movie, Colormap] = eegmovie(data,srate,eloc_locs,titl,movieframes,minmax,startsec,varargin)

if nargin<1
	help eegmovie
    return
end

clf
[chans,frames] = size(data);

icadefs;   % read DEFAULT_SRATE;

if nargin<7
   startsec = 0;
end
if nargin<6
   minmax = 0;
end
if minmax ==0,
	datamin = min(min(data));
	datamax = max(max(data));
    absmax  = max([abs(datamin), abs(datamax)]);
    fudge   = 0.05*(datamax-datamin); % allow for slight extrapolation
    datamin = -absmax-fudge;
    datamax =  absmax+fudge;
    minmax = [datamin datamax];
end
if nargin <5
	movieframes = 0;
end
if movieframes == 0
	movieframes = 1:frames;
end
if nargin <4
	titl = '';
end
if titl == 0
	titl = '';
end
if nargin <2
	srate = 0;
end
if nargin <3
	eloc_locs = 0;
end

if movieframes(1) < 1 | movieframes(length(movieframes))>frames
	fprintf('eegmovie(): specified movieframes not in data!\n');
	return
end

if srate ==0,
	srate = DEFAULT_SRATE;
end

mframes = length(movieframes);
fprintf('Making a movie of %d frames\n',mframes)
Movie    = moviein(mframes,gcf);

%%%%%%%%%%%%%%%%%%%%% eegplot() of data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axeegplot = axes('Units','Normalized','Position',[.75 .05 .2 .9]);

% >> eegplotold('noui',data,srate,spacing,eloc_file,startsec,color)
if isstruct(eloc_locs)
    eegplotold('noui',-data,srate,0,[],startsec,'r');
else
    eegplotold('noui',-data,srate,0,eloc_locs,startsec,'r');
end;    
% set(axeegplot,'XTick',[])                %%CJH
                                         % plot negative up
limits = get(axeegplot,'Ylim');          % list channel numbers only
set(axeegplot,'GridLineStyle',':')
set(axeegplot,'Xgrid','off')
set(axeegplot,'Ygrid','on')
axcolor = get(gcf,'Color');

%%%%%%%%%%%%%%%%%%%%% topoplot() maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axtopoplot = axes('Units','Normalized','Position',[0 .1 .72 .8],'Color',axcolor);
TPAXCOLOR  = get(axtopoplot,'Color');    %%CJH
Colormap   = [jet(64);TPAXCOLOR];        %%CJH

fprintf('Warning: do not resize plot window during movie creation ...\n   ');

%%%%%%%%%%%%%%%%%%%%%%%%% "Roll'em!" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for f = 1:length(movieframes)                      % make the movie, frame by frame
   i=movieframes(f);
   fprintf('%d ',f);

   axes(axeegplot)
   x1 = startsec+(i-1)/srate;
   l1 = line([i i],limits,'color','b'); % draw vertical line at map timepoint
   timetext = num2str(x1,3);   % Note: Matlab4 doesn't take '%4.3f'
   set(axeegplot,'Xtick',i,'XtickLabel',num2str(x1,'%4.3f'));

   axes(axtopoplot)
   cla
   set(axtopoplot,'Color',axcolor);
   topoplot(data(:,i),eloc_locs,'style','both','maplimits',minmax, varargin{:}); 
					                     % use channel locations file
   txt = [ int2str(f)];  
   text(-0.5,-0.5,txt,'FontSize',14);    % show frame number
   title(titl,'FontSize',16)

   Movie(:,f) = getframe(gcf);
   drawnow
   delete(l1)

   if rem(i,8) == 0
     fprintf('\n'); % print newlines 
   end
end
fprintf('\nDone\n');
