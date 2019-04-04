% compmap() - Plot multiple topoplot() maps of ICA component topographies
%             Click on an individual map to view separately. 
% Usage:
%       >> compmap (winv,'eloc_file',compnos,'title',rowscols,labels,printflag)
%
% Inputs:
%   winv       - Inverse weight matrix = EEG scalp maps. Each column is a
%                map; the rows correspond to the electrode positions
%                defined in the eloc_file. Normally, winv = inv(weights*sphere).
%  'eloc_file' - Name of the eloctrode position file in the style defined
%                by >> topoplot example {default|0 ->'chan_file'}
%   compnos    - Vector telling which (order of) component maps to show
%                Indices <0 tell compmap to invert a map; = 0 leave blank sbplot 
%                Example [1 0 -2 3 0 -6] {default|0 -> 1:columns_in_winv}
%  'title'     - Title string for each page {default|0 -> 'ICA Component Maps'}
%   rowscols   - Vector of the form [m,n] where m is total vertical tiles and n 
%                is horizontal tiles per page. If the number of maps exceeds m*n,
%                multiple figures will be produced {def|0 -> one near-square page}.
%   labels     - Vector of numbers or a matrix of strings to use as labels for
%                each map, else ' ' -> no labels {default|0 -> 1:ncolumns_in_winv}
%   printflag  - 0= screen-plot colors {default}
%                1= printer-plot colors
%
% Note: Map scaling is to +/-max(abs(data); green = 0
%
% Author: Colin Humphries, CNL / Salk Institute, Aug, 1996
%
% See also: topoplot()

% This function calls topoplot(). and cbar().

% Copyright (C) Colin Humphries & Scott Makeig, CNL / Salk Institute, Aug, 1996
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

% Colin Humphries, CNL / Salk Institute, Aug, 1996
%    03-97 revised -CH
% 11-05-97 revised for Matlab 5.0.0.4064; added negative compnnos option
%          improved help msg; added figure offsets -Scott Makeig & CH
% 11-13-97 fixed compnos<0 bug -sm
% 11-18-97 added test for labels=comps -sm 
% 12-08-97 added test for colorbar_tp() -sm
% 12-15-97 added axis('square'), see SQUARE below -sm
% 03-09-98 added test for eloc_file, fixed size checking for labels -sm
% 02-09-00 added test for ' ' for srclabels(1,1) -sm
% 02-23-00 added srclabels==' ' -> no labels -sm
% 03-16-00 added axcopy() -sm & tpj
% 02-25-01 added test for max(compnos) -sm
% 05-24-01 added HEADPLOT logical flag below -sm
% 01-25-02 reformated help & license -ad 

% NOTE: 
% There is a minor problem with the Matlab function colorbar().
% Use the toolbox cbar() instead.

function compmap(Winv,eloc_file,compnos,titleval,pagesize,srclabels,printlabel,caxis)


DEFAULT_TITLE = 'ICA Component Maps';
DEFAULT_EFILE = 'chan_file';
NUMCONTOUR = 5;     % topoplot() style settings
OUTPUT = 'screen';  % default: 'screen' for screen colors, 
                    %          'printer' for printer colors
STYLE = 'both';
INTERPLIMITS = 'head';
MAPLIMITS    = 'absmax';
SQUARE       = 1; % 1/0 flag making topoplot() asex square -> round heads
ELECTRODES   = 'on'; % default: 'on' or 'off'
ELECTRODESIZE  = []; % defaults 1-10 set in topoplot text.
HEADPLOT     = 0; % 1/0 plot 3-D headplots instead of 2-d topoplots.

if nargin<1
   help compmap
   return
end

curaxes = gca;
curpos = get(curaxes,'Position');
mapaxes = axes('Position',[curpos(1) curpos(2)+0.09 curpos(3) curpos(4)-0.09],...
                     'Visible','off'); 
                             % leave room for cbar
set(mapaxes,'visible','off');
pos = get(mapaxes,'Position');

% delete(gca);
% ax = axes('Position',pos);

[chans, frames] = size (Winv);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check inputs and set defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==8
  if strcmp(caxis(1,1:2), 'mi') % min/max of data
     MAPLIMITS = [min(min(Winv(:,compnos))) max(max(Winv(:,compnos)))];
  elseif caxis(1,1:2) == 'ab' % +/-max abs data
     absmax = max([abs(min(min(Winv(:,compnos)))) abs(max(max(Winv(:,compnos))))]);
     MAPLIMITS = [-absmax absmax];
  elseif size(caxis) == [1,2] % given color axis limits
     MAPLIMITS = caxis;
  end % else default
end
if nargin < 7
   printlabel = 0;
end
if printlabel == 0
   printlabel = OUTPUT; % default set above
else
   printlabel = 'printer';
end

if nargin < 6
   srclabels = 0;
end
if nargin < 5
   pagesize = 0;  
end
if nargin < 4
   titleval = 0;
end
if nargin < 3
   compnos = 0;
end
if nargin < 2
   eloc_file = 0;
end

if srclabels == 0
  srclabels = [];
end
if titleval == 0;
   titleval = DEFAULT_TITLE;
end
if compnos == 0
   compnos = (1:frames);
end
if max(compnos)>frames
   fprintf('compmap(): Cannot show comp %d. Only %d components in inverse weights\n',...
      max(compnos),frames);
   return
end
if pagesize == 0
   numsources = length(compnos);
   DEFAULT_PAGE_SIZE = ...
[floor(sqrt(numsources)) ceil(numsources/floor(sqrt(numsources)))];
   m = DEFAULT_PAGE_SIZE(1);
   n = DEFAULT_PAGE_SIZE(2);
elseif length(pagesize) ==1
   help compmap
   return
else
   m = pagesize(1);
   n = pagesize(2);
end
if eloc_file == 0 
   eloc_file = DEFAULT_EFILE;
end

totalsources = length(compnos);
if ~isempty(srclabels) 
  if ~ischar(srclabels(1,1)) || srclabels(1,1)==' ' % if numbers
    if size(srclabels,1) == 1
       srclabels = srclabels';
    end
  end
  if size(srclabels,1)==1 && size(srclabels,2)==1 && srclabels==' ' 
     srclabels = repmat(srclabels,totalsources,1);
  end
  if size(srclabels,1) ~= totalsources,
     fprintf('compmap(): numbers of components and component labels do not agree.\n');
     return
  end
end
pages = ceil(totalsources/(m*n));		
if pages > 1
   fprintf('compmap(): will create %d figures of %d by %d maps: ',...
            pages,m,n);
end

off = [ 25 -25 0 0];  % position offsets for multiple figures

fid = fopen(eloc_file);
if fid<1,
  fprintf('compmap()^G: cannot open eloc_file (%s).\n',eloc_file);
  return
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot the maps %%%%%%%%%%%%%%%%%%%%%%%

for i = (1:pages)
  if i > 1
    figure('Position',pos+(i-1)*off); % place figures in right-downward stack
    set(gca,'Color','w') %CJH - set background color to white
    curaxes = gca;
    curpos = get(curaxes,'Position'); % new whole-figure axes
  end
  
  if (totalsources > i*m*n)
    sbreak = n*m;
  else 
    sbreak = totalsources - (i-1)*m*n; % change page/figure after this many 
  end

  for j = (1:sbreak) % maps on this page/figure
    comp = j+(i-1)*m*n; % compno index
    if compnos(comp)~=0
      if compnos(comp)>0
       source_var = Winv(:,compnos(comp))';       % plot map
      elseif compnos(comp)<0
       source_var = -1*Winv(:,-1*compnos(comp))'; % invert map
      end

      sbplot(m,n,j,'ax',mapaxes);
      if HEADPLOT
       headplot(source_var,eloc_file,'electrodes','off'); % 3-d image
      else
       topoplot(source_var,eloc_file,...
        'style',STYLE,...
        'electrodes',ELECTRODES,...
        'emarkersize',ELECTRODESIZE,...
        'numcontour',NUMCONTOUR,...
        'interplimits',INTERPLIMITS,...
        'maplimits',MAPLIMITS);             % draw 2-d scalp map
      end
      if SQUARE,
         axis('square');
      end

      if isempty(srclabels)
        title(int2str(compnos(comp)))   ;
      else
        if ischar(srclabels)      
          title(srclabels(comp,:));
        else
	      title(num2str(srclabels(comp)));
        end
      end
      drawnow  % draw one map at a time
    end
  end

  % ax = axes('Units','Normal','Position',[.5 .06 .32 .05],'Visible','Off');
  axes(curaxes);
  set(gca,'Visible','off','Units','normalized');
  curpos = get(gca,'position');
  ax = axes('Units','Normalized','Position',...
      [curpos(1)+0.5*curpos(3) curpos(2)+0.01*curpos(4) ...
       0.32*curpos(3) 0.05*curpos(4)],'Visible','Off');
  if exist('cbar') == 2
    cbar(ax);        % Slightly altered Matlab colorbar()
                     % Write authors for further information.
  else
    colorbar(ax);    % Note: there is a minor problem with this call.
  end
  axval = axis;
  Xlim = get(ax,'Xlim');
  set(ax,'XTick',(Xlim(2)+Xlim(1))/2);
  set(ax,'XTickMode','manual');
  set(ax,'XTickLabelMode','manual');
  set(ax,'XTickLabel','0');

  axes(curaxes);
  axis off;
  % axbig = axes('Units','Normalized','Position',[0 0 1 1],'Visible','off');
  t1 = text(.25,.07,titleval,'HorizontalAlignment','center');
  if pages > 1
     fprintf('%d ',i);
  end
  axcopy(gcf); % allow popup window of single map with mouse click
end
if pages > 1
   fprintf('\n');
end
