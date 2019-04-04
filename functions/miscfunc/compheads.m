% compheads() - plot multiple topoplot() maps of ICA component topographies
%
% Usage:
%       >> compheads(winv,'spline_file',compnos,'title',rowscols,labels,view)
%
% Inputs:
%   winv       - Inverse weight matrix = EEG scalp maps. Each column is a
%                map; the rows correspond to the electrode positions
%                defined in the eloc_file. Normally, winv = inv(weights*sphere).
%  spline_file - Name of the eloctrode position file in BESA spherical coords.
%   compnos    - Vector telling which (order of) component maps to show
%                Indices <0 tell compheads() to invert a map; = 0 leave blank subplot 
%                Example [1 0 -2 3 0 -6] {default|0 -> 1:columns_in_winv}
%  'title'     - Title string for each page {default|0 -> 'ICA Component Maps'}
%   rowscols   - Vector of the form [m,n] where m is total vertical tiles and n 
%                is horizontal tiles per page. If the number of maps exceeds m*n,
%                multiple figures will be produced {def|0 -> one near-square page}.
%   labels     - Vector of numbers or a matrix of strings to use as labels for
%                each map {default|0 -> 1:ncolumns_in_winv}
%   view       - topoplot() view, either [az el] or keyword ('top',...)
%                See >> help topoplot() for options.
%
% Note: Map scaling is to +/-max(abs(data); green = 0
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 4-28-1998 
%
% See also: topoplot()

% Copyright (C) 4-28-98 from compmap.m Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

function compheads(Winv,eloc_file,compnos,titleval,pagesize,srclabels,view)

if nargin<1
   help compheads
   return
end

[chans, frames] = size (Winv);

DEFAULT_TITLE = 'ICA Component Maps';
DEFAULT_EFILE = 'chan_file';
NUMCONTOUR = 5;     % topoplot() style settings
OUTPUT = 'screen';  % default: 'screen' for screen colors, 
                    %          'printer' for printer colors
STYLE = 'both';
INTERPLIMITS = 'head';
MAPLIMITS = 'absmax';
SQUARE    = 1; % 1/0 flag making topoplot() asex square -> round heads

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check inputs and set defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

printlabel = OUTPUT; % default set above
if nargin < 7
   view = [-127 30];
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
if pagesize == 0
   numsources = length(compnos);
   DEFAULT_PAGE_SIZE = ...
[floor(sqrt(numsources)) ceil(numsources/floor(sqrt(numsources)))];
   m = DEFAULT_PAGE_SIZE(1);
   n = DEFAULT_PAGE_SIZE(2);
elseif length(pagesize) ==1
   help compheads
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
  if ~ischar(srclabels(1,1)) % if numbers
    if size(srclabels,1) == 1
       srclabels = srclabels';
    end
  end
  if size(srclabels,1) ~= totalsources,
     fprintf('compheads(): numbers of components and component labels do not agree.\n');
     return
  end
end
pages = ceil(totalsources/(m*n));		
if pages > 1
   fprintf('compheads(): will create %d figures of %d by %d maps: ',...
            pages,m,n);
end

pos = get(gcf,'Position');
off = [ 25 -25 0 0];  % position offsets for multiple figures

fid = fopen(eloc_file);
if fid<1,
  fprintf('compheads()^G: cannot open eloc_file (%s).\n',eloc_file);
  return
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot the maps %%%%%%%%%%%%%%%%%%%%%%%

for i = (1:pages)
  if i > 1
    figure('Position',pos+(i-1)*off); % place figures in right-downward stack
  end
  set(gcf,'Color','w') %CJH - set background color to white
  
  if (totalsources > i*m*n)
    sbreak = n*m;
  else 
    sbreak = totalsources - (i-1)*m*n;
  end

  for j = (1:sbreak) % maps on this page
    comp = j+(i-1)*m*n; % compno index
    if compnos(comp)~=0
      if compnos(comp)>0
       source_var = Winv(:,compnos(comp))';       % plot map
      elseif compnos(comp)<0
       source_var = -1*Winv(:,-1*compnos(comp))'; % invert map
      end

      subplot(m,n,j)
      headplot(source_var,eloc_file,'electrodes','off','view',view);
      %topoplot(source_var,eloc_file,'style',STYLE,...
      % 'numcontour',NUMCONTOUR,'interplimits',INTERPLIMITS,...
      % 'maplimits',MAPLIMITS); % draw map
      if SQUARE,
         axis('square');
      end

      if isempty(srclabels)
          t=title(int2str(compnos(comp)));
          set(t,'FontSize',16);
      else
        if ischar(srclabels)      
          t=title(srclabels(comp,:));
          set(t,'FontSize',16);
        else
	      t=title(num2str(srclabels(comp)));
          set(t,'FontSize',16);
        end
      end
      drawnow  % draw one map at a time
    end
  end

  ax = axes('Units','Normal','Position',[.5 .04 .32 .05],'Visible','Off');
  colorbar(ax)     
  axval = axis;
  Xlim = get(ax,'Xlim');
  set(ax,'XTick',(Xlim(2)+Xlim(1))/2)
  set(gca,'XTickLabel','0')
  set(gca,'XTickMode','manual')

  axbig = axes('Units','Normalized','Position',[0 0 1 1],'Visible','off');
  t1 = text(.25,.070,titleval,'HorizontalAlignment','center','FontSize',14);
  if pages > 1
     fprintf('%d ',i);
  end
end
if pages > 1
   fprintf('\n');
end
