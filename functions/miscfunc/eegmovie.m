% eegmovie() - Compile and view a Matlab movie. 
%              Uses scripts eegplotold() and topoplot().
%              Use seemovie() to display the movie.
% Usage:
% >> [Movie,Colormap] = eegmovie(data,srate,elec_locs, 'key', val, ...);
%
% Or legacy call
% >> [Movie,Colormap] = eegmovie(data,srate,elec_locs,title,movieframes,minmax,startsec,...);
%
% Inputs:
%   data        = (chans,frames) EEG data set to plot
%   srate       = sampling rate in Hz
%   elec_locs   = electrode locations structure or file
%
% Optional inputs:
%   'mode'        = ['2D'|'3D'] plot in 2D using topoplot or in 3D using
%                   headplot. Default is 2D.
%   'headplotopt' = [cell] optional inputs for headplot. Default is none.
%   'topoplotopt' = [cell] optional inputs for topoplot. Default is none.
%   'title'       = plot title. Default is none.
%   'movieframes' = vector of frames indices to animate. Default is all.
%   'minmax'      = [blue_lower_bound, red_upper_bound]. Default is 
%                   +/-abs max of data.
%   'startsec'    = starting time in seconds. Default is 0.
%   'timecourse'  = ['on'|'off'] show time course for all electrodes. Default is 'on'.
%   'framenum'    = ['on'|'off'] show frame number. Default is 'on'.
%   'time'        = ['on'|'off'] show time in ms. Default is 'off'.
%   'vert'        = [float] plot vertical lines at given latencies. Default is none.
%   'camerapath'  = [az_start az_step el_start el_step] {default [-127 0 30 0]}
%                   Setting all four non-0 creates a spiral camera path
%                   Subsequent rows [movieframe az_step 0 el_step] adapt step 
%                   sizes allowing starts/stops, panning back and forth, etc.
%
% Legacy inputs:
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
% Author: Arnaud Delorme, Colin Humphries & Scott Makeig, CNL, Salk Institute, La Jolla, 3/97
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

function [Movie, Colormap] = eegmovie(data,srate,eloc_locs,varargin);
%titl,movieframes,minmax,startsec,varargin)

if nargin<1
	help eegmovie
    return
end

clf
[chans,frames] = size(data);
icadefs;   % read DEFAULT_SRATE;
if nargin <2
	srate = 0;
end
if nargin <3
	eloc_locs = 0;
end

if nargin > 5 && ~ischar(varargin{3}) || nargin == 4 && ~ischar(varargin{2})
    % legacy mode
    options = {};
    if nargin>=8, options = { options{:} 'topoplotopt' varargin(5:end) }; end
    if nargin>=7, options = { options{:} 'startsec'    varargin{4} }; end
    if nargin>=6, options = { options{:} 'minmax'      varargin{3} }; end
    if nargin>=5, options = { options{:} 'movieframes' varargin{2} }; end
    if nargin>=4, options = { options{:} 'title'       varargin{1} }; end
else
    options = varargin;
end;

opt = finputcheck(options, { 'startsec'    'real'    {}    0;
                             'minmax'      'real'    {}    0;
                             'movieframes' 'integer' {}    0;
                             'title'       'string'  {}    '';
                             'vert'        'real'    {}    [];
                             'mode'        'string'  { '2D' '3D'  }    '2D';
                             'timecourse'  'string'  { 'on' 'off' }    'on';
                             'framenum'    'string'  { 'on' 'off' }    'on';
                             'camerapath'  'real'    []                0;
                             'time'        'string'  { 'on' 'off' }    'off';
                             'topoplotopt' 'cell'    {}    {};
                             'headplotopt' 'cell'    {}    {} }, 'eegmovie');
if isstr(opt), error(opt); end;
if opt.minmax ==0,
	datamin = min(min(data));
	datamax = max(max(data));
    absmax  = max([abs(datamin), abs(datamax)]);
    fudge   = 0.05*(datamax-datamin); % allow for slight extrapolation
    datamin = -absmax-fudge;
    datamax =  absmax+fudge;
    opt.minmax = [datamin datamax];
end

if opt.movieframes == 0
	opt.movieframes = 1:frames;
end
if opt.movieframes(1) < 1 || opt.movieframes(length(opt.movieframes))>frames
	fprintf('eegmovie(): specified movieframes not in data!\n');
	return
end
if srate ==0,
	srate = DEFAULT_SRATE;
end
if strcmpi(opt.time, 'on'), opt.framenum = 'off'; end;

mframes = length(opt.movieframes);
fprintf('Making a movie of %d frames\n',mframes)
Movie    = moviein(mframes,gcf);

%%%%%%%%%%%%%%%%%%%%% eegplot() of data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(opt.timecourse, 'on')
    axeegplot = axes('Units','Normalized','Position',[.75 .05 .2 .9]);
    
    % >> eegplotold('noui',data,srate,spacing,eloc_file,startsec,color)
    if isstruct(eloc_locs)
        fid = fopen('tmp_file.loc', 'w');
        adddots = '...';
        for iChan = 1:length(eloc_locs)
            fprintf(fid, '0 0 0 %s\n', [ eloc_locs(iChan).labels adddots(length(eloc_locs(iChan).labels):end) ]);
        end;
        fclose(fid);
        eegplotold('noui',-data,srate,0,'tmp_file.loc',opt.startsec,'r');
    else
        eegplotold('noui',-data,srate,0,eloc_locs,opt.startsec,'r');
    end;
    
    % set(axeegplot,'XTick',[])                %%CJH
    % plot negative up
    limits = get(axeegplot,'Ylim');          % list channel numbers only
    set(axeegplot,'GridLineStyle',':')
    set(axeegplot,'Xgrid','off')
    set(axeegplot,'Ygrid','on')
    
    for ind = 1:length(opt.vert)
       frameind = (opt.vert(ind)-opt.startsec)*srate+1;
       line([frameind frameind],limits,'color','k'); % draw vertical line at map timepoint
       set(axeegplot,'Xtick',frameind,'XtickLabel',num2str(opt.vert(ind),'%4.3f'));
    end;
end;

%%%%%%%%%%%%%%%%%%%%% topoplot/headplot axis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axcolor = get(gcf,'Color');
axtopoplot = axes('Units','Normalized','Position',[0 .1 .72 .8],'Color',axcolor);
TPAXCOLOR  = get(axtopoplot,'Color');    %%CJH
Colormap   = [jet(64);TPAXCOLOR];        %%CJH
fprintf('Warning: do not resize plot window during movie creation ...\n   ');
h = textsc(opt.title, 'title'); set(h,'FontSize',16)

%%%%%%%%%%%%%%%%%%%%% headplot setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(opt.mode, '3d')
    headplot('setup',eloc_locs, 'tmp.spl', opt.headplotopt{:});
    
    if isequal(opt.camerapath, 0)
        opt.camerapath = [-127 0 30 0];
        fprintf('Using default view [-127 0 30 0].');
    end;
    if size(opt.camerapath,2)~=4
        error('Camerapath parameter must have exact 4 columns');
    end

    newpan = length(opt.movieframes)+1;
    posrow = 2;
    if size(opt.camerapath,1) > 1
        newpan = opt.camerapath(posrow,1); % pick up next frame to change camerapos step values
    end
    
    azimuth   = opt.camerapath(1,1); % initial camerapath variables
    az_step   = opt.camerapath(1,2);
    elevation = opt.camerapath(1,3);
    el_step   = opt.camerapath(1,4);
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%% "Roll'em!" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for f = 1:length(opt.movieframes)                      % make the movie, frame by frame
   indFrame = opt.movieframes(f);

   % show time course
   if strcmpi(opt.timecourse, 'on')
       axes(axeegplot)
       x1 = opt.startsec+(indFrame-1)/srate;
       l1 = line([indFrame indFrame],limits,'color','b'); % draw vertical line at map timepoint
       set(axeegplot,'Xtick',indFrame,'XtickLabel',num2str(x1,'%4.3f'));
   end;
   
   % plot headplot or topoplot
   axes(axtopoplot)
   cla
   set(axtopoplot,'Color',axcolor);
   if strcmpi(opt.mode, '2d')
       topoplot(data(:,indFrame),eloc_locs,'maplimits',opt.minmax, opt.topoplotopt{:}); 
   else
       headplot(data(:,indFrame),'tmp.spl','view',[azimuth elevation], opt.headplotopt{:});
       
       % adapt camerapath step sizes
       if indFrame == newpan   
           az_step = opt.camerapath(posrow,2);
           el_step = opt.camerapath(posrow,4);
           posrow = posrow+1;
           if size(opt.camerapath,1)>=posrow
               newpan = opt.camerapath(posrow,1);
           else
               newpan = length(opt.movieframes)+1;
           end
       end
       
       % update camera position
       azimuth = azimuth+az_step;     
       elevation = elevation+el_step;
       if elevation>=90
           fprintf('headplot(): warning -- elevation out of range!\n');
           elevation = 89.99;
       end
       if elevation<=-90
           fprintf('headplot(): warning -- elevation out of range!\n');
           elevation = -89.99;
       end
       
       set(axtopoplot,'Units','pixels',...
           'CameraViewAngleMode','manual',...
           'YTickMode','manual','ZTickMode','manual',...
           'PlotBoxAspectRatioMode','manual',...
           'DataAspectRatioMode','manual');    % keep camera distance constant
       
   end;
                                         
   % show frame number
   if strcmpi(opt.framenum, 'on') 
       txt = [ int2str(f)]; 
       text(-0.5,-0.5,txt,'FontSize',14);    
   elseif strcmpi(opt.time, 'on') 
       txt = sprintf('%3.3f s', opt.startsec+(indFrame-1)/srate); 
       text(-0.5,-0.5,txt,'FontSize',14);    
   end;

   Movie(:,f) = getframe(gcf);
   drawnow
   if strcmpi(opt.timecourse, 'on')
       delete(l1)
   end;
   
   % print advancement
   fprintf('.',f);
   if rem(indFrame,10) == 0, fprintf('%d',f); end;
   if rem(indFrame,50) == 0, fprintf('\n'); end
end
fprintf('\nDone\n');
