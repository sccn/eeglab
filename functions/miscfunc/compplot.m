% compplot() - plot a data epoch and maps its scalp topography at a given time
%
% Usage: To plot the projection of an ICA component onto the scalp
%        >> projdata = icaproj(data,weights,compindex);
%
% then   >> compplot(projdata);
%
% else to plot an EEG epoch with a topoplot at one selected time point
%        >> compplot(data,plotframe,chan_file,xstart,srate,title, splinefile);
%
% Inputs:
%  data        = data returned from icaproj() *ELSE* any EEG/ERP data epoch
%  plotframe   = frame to plot topographically {default|0 -> frame of max(var())}
%  'chan_file' = chan file, see >> topoplot example {def|0 -> 'chan_file'}
%  xstart      = start time in seconds {default|0 -> 0}
%  srate       = data sampling rate in Hz {default|0 -> 256 Hz}
%  'title'     = plot title {default|0 -> none}
%  splinefile  = headplot spline file (optional) for 3-d head image {default: none}
%
% Authors: Colin Humphries & Scott Makeig, SCCN/INC/UCSD, CNL / Salk Institute, 1997 
%
% See also: icaproj()

% Copyright (C) 2000 Colin Humphries & Scott Makeig, SCCN/INC/UCSD
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

% 6-07-97 changed order of args in icaproj above -sm
% Revised for Matlab 5.0.0.4064 11/97 -ch & sm
% 2-18-98 added undocumented headplot() option -sm
% 3-09-98 added check for chan_file -sm
% 3-18-98 use new eegplot('noui') -sm
% 7-25-98 make sure length(hkids) >= 3 on line 122 ff -sm
% 12-19-00 updated icaproj() call in help msg -sm
% 1-24-01 corrected input error in header -ad
% 01-25-02 reformated help & license, added links -ad 

function M = compplot(data,plotframe,chan_file,xstart,srate,titl,spline_file)

if nargin < 1
   help compplot
   return
end

[chans,frames] = size(data);
icadefs;   % read DEFAULT_SRATE

if nargin < 7
    spline_file = nan;
end
if nargin < 6,
   titl = '';     % DEFAULT TITLE
end
if nargin < 5
    srate = 0;  
end
if nargin < 4
    xstart = 0;   % DEFAULT START TIME
end
if nargin < 3
    chan_file = 'chan_file';  % DEFAULT CHAN_FILE
end
if chan_file == 0,
    chan_file = 'chan_file';  % DEFAULT CHAN_FILE
end
if nargin < 2
    plotframe = 0;
end
if plotframe == 0
    [mx plotframe] = max(mean(data.*data)); 
                         % default plotting frame has max variance
end
if plotframe > frames
    fprintf('Plot frame %d is > frames in data (%d)\n',plotframe,frames);
    return
end
fprintf('Topoplot will show frame %d\n',plotframe);
if nargin < 1
    help compplot
end

if srate==0,
	srate = DEFAULT_SRATE;
end

clf % clear the current figure
set(gcf,'Color',BACKCOLOR); % set the background color
axe = axes('Units','Normalized','Position',[.75 .12 .2 .8]);
set(gca,'Color','none');

% >> eegplot('noui',data,srate,spacing,eloc_file,startsec,color)
eegplot('noui',-data(:,xstart*srate+1:end),srate,0,chan_file,0,'r')

plottime = xstart + (plotframe-1)/srate;
timetext = num2str(plottime,'%4.3f');
set(gca,'XTick',plotframe,'XtickLabel',timetext)       %%CJH
set(gca,'Color','none');

limits = get(gca,'Ylim');
set(gca,'GridLineStyle',':')
set(gca,'Xgrid','off')
set(gca,'Ygrid','on')

% axes(axe)
plottime = xstart + (plotframe-1)/srate;
l1 = line([plotframe plotframe],limits,'color','b'); % was [plottime plottime]

fid = fopen(chan_file); % check whether topoplot will find chan_fil
if fid<1,
  fprintf('topoplot()^G: cannot open eloc_file (%s).\n',chan_file);
  return
end

axcolor = get(gcf,'Color');
axt = axes('Units','Normalized','Position',[.00 .10 .75 .80],'Color',axcolor); 
axes(axt)                                                       % topoplot axes
cla
if isnan(spline_file)
   topoplot(data(:,plotframe),chan_file,'style','both');
else
   headplot(data(:,plotframe),spline_file);
end
text(0.00,0.70,titl,'FontSize',14,'HorizontalAlignment','Center');

axt = axes('Units','Normalized','Position',[.05 .05 .055 .18]); % colorbar axes
h=cbar(axt);  

% hkids=get(h,'children');
% delno = 3;
% if delno > length(hkids)
%    delno = length(hkids);
% end
% for i=1:delno
%    delete(hkids(i));
% end

axis(axis)
set(h,'XTickLabel',[]);
HYLIM = get(h,'YLim');              %%CJH
set(h,'YTick',[HYLIM(1) (HYLIM(1)+HYLIM(2))/2 HYLIM(2)])  %%CJH
set(h,'YTickLabel',['-';'0';'+']);     %%CJH
set(h,'YTickMode','manual');

axt = axes('Units','Normalized','Position',[0 0 1 1],...
            'Visible','Off','Fontsize',16);                     % topoplot axes
axes(axt)
timetext = [ 't = ',num2str(plottime) ' s'];
text(0.16,0.09,timetext,'FontSize',14);
text(0.85,0.04,'Time (sec)','FontSize',14,'HorizontalAlignment','Center');

