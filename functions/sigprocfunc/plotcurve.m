% plotcurve() - plot curve(s) with optional significance highlighting.
%
% Usage: >> plotcurve(times, data, 'key1', 'val1', 'key2', val2' ...);
%
% Required inputs:
%   times - [float] vector of time indices
%   data  - [float] data array, size of [n x times]. If n>1 several
%           curves are plotted (unless 'plotmean' option is used).
%
% Optional inputs:
%  'maskarray' = Input bootstrap limits. Can be 1-D [min max], 2-D (min
%                and max for all ordinate) or 3-D (min and max for all
%                ordinates at all time points).
%  'val2mask'  = Value to use for generating mask. By default use data.
%  'plotmean' = ['on'|'off'] For 'curve' plots only. Average all
%                frequencies given as input. Default: 'on'.
%  'highlightmode'  = ['background'|'bottom'] For 'curve' plots only,
%                display significant time regions either in the plot
%                background or underneatht the curve.
%  'xlabel'    = [string] x label
%  'ylabel'    = [string] y label
%  'legend'    = [cell] legend. Cell array of string, with one string
%                per curve.
%  'ylim'      = [min max] or [min] limits for the ordinate axis.
%  'title'     = Optional figure title. If two conditions are given
%                as input, title can be a cell array with two text
%                string elements {none}
%  'vert'      = Latencies to mark with a dotted vertical line   {none}
%  'linewidth' = Line width for marktimes traces (thick=2, thin=1) {2}
%  'chanlocs'  = channel location structure.
%  'plottopo'  = [min max] plot topography within the time limits defined
%                in this function. If several lines are given as input, one
%                scalp map is plot for each line.
%  'traceinfo' = [string|cell array] information shown on the command line
%                when the user click on a specific trace. Default none.
%
% Authors: Arnaud Delorme, 2004, Bhaktivedanta Institute

% Copyright (C) 2004  Arnaud Delorme
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

function plotcurve( times, R, varargin);

	if nargin < 2
        help plotcurve;
        return;
	end
	g = finputcheck( varargin, { 'maskarray' ''        []       [];
                             'val2mask'      'real'    []           R;
                             'highlightmode' 'string'  { 'background','bottom' } 'background';
                             'plotmean'      'string'  { 'on','off' }            'off';
                             'plotindiv'     'string'  { 'on','off' }            'on';
                             'traceinfo'   { 'string' 'cell' } { { } {} }        'off';
                             'logpval'       'string'  { 'on','off' }            'off';
                             'title'         'string'  []                        '';
                             'xlabel'        'string'  []                        '';
                             'plotmode'      'string'  {'single','topo'}         'single';
                             'ylabel'        'string'  []                        '';
                             'legend'        'cell'    []                        {};
                             'colors'        'cell'    []                        {};
                             'plottopotitle' 'cell'    []                        {};
                             'chanlocs'      'struct'  []                        struct;
                             'ylim'          'real'    []                        [];
                             'vert'          'real'    []                        [];
                             'plottopo'      'real'    []                        [];
                             'linewidth'     'real'    []                        2;
                             'marktimes'     'real'    []                        [] });
   if ischar(g), error(g); end
  % keyboard;
   if isempty(g.colors), g.colors = { 'r' 'g' 'b' 'c' 'm' 'r' 'b' 'g' 'c' 'm' 'r' 'b' 'g' 'c' 'm' 'r' 'b' ...
                   'g' 'c' 'm' 'r' 'b' 'g' 'c' 'm' 'r' 'b' 'g' 'c' 'm' 'r' 'b' 'g' 'c' 'm' }; end
   if strcmpi(g.plotindiv, 'off'), g.plotmean = 'on'; end
   
   if ~any(length(times) == size(R))
       try
           R = reshape(R, length(times), length(R)/length(times))';
       catch, error('Size of time input and array input does not match');
       end
   else
       if size(R,3) > 1 && size(R,1) == 1
           R = permute(R, [3 2 1]);
       end
   end
   
   % regions of significance
   % -----------------------
   if ~isempty(g.maskarray)
       if length(unique(g.maskarray)) < 4
          Rregions = g.maskarray;
       else
          Rregions = ones(size(g.val2mask));
          switch dims(g.maskarray)
               case 3, Rregions  (find(g.val2mask > g.maskarray(:,:,1) & (g.val2mask < g.maskarray(:,:,2)))) = 0;
               case 2, if size(g.val2mask,2) == size(g.maskarray,2)
                           Rregions  (find(g.val2mask < g.maskarray)) = 0;
                       elseif size(g.val2mask,1) == size(g.maskarray,1)
                           Rregions(find((g.val2mask > repmat(g.maskarray(:,1),[1 length(times)])) ...
                             & (g.val2mask < repmat(g.maskarray(:,2),[1 length(times)])))) = 0;
                       else
                           Rregions(find((g.val2mask > repmat(g.maskarray(:,1),[1 length(times)])) ...
                             & (g.val2mask < repmat(g.maskarray(:,2),[1 length(times)])))) = 0;
                       end
               case 1, Rregions  (find(g.val2mask < repmat(g.maskarray(:),[1 size(g.val2mask,2)]))) = 0;
           end; 
           Rregions = sum(Rregions,1);
       end
   else 
       Rregions = [];
   end

  % plotting
  % --------
  if size(R,1) == length(times), R = R'; end
  if strcmpi(g.plotmean, 'on') || strcmpi(g.plotindiv, 'off')
      if strcmpi(g.plotindiv, 'on')
          R = [ R; mean(R,1) ];
      else
          R = mean(R,1);
      end
  end
  ax = gca;
  if ~isempty(g.maskarray) && strcmpi(g.highlightmode, 'bottom')
      pos = get(gca, 'position');
      set(gca, 'position', [ pos(1)+pos(3)*0.1 pos(2)+pos(4)*0.1 pos(3)*0.9 pos(4)*0.85 ]);
  end
  
  % plot topographies
  % -----------------
  if ~isempty(g.plottopo)
      tmpax = gca;
      pos = get(gca, 'position');
      set(gca, 'position', [ pos(1) pos(2) pos(3) pos(4)/2 ]);
      
      for index = 1:size(g.plottopo)
          axes('position', [ (index-1)*pos(3)/size(g.plottopo,1)+pos(1) pos(2)+pos(4)/2 pos(3)/size(g.plottopo,1) pos(4)/2 ]);
          %topoplot(g.plottopo(index,:), g.chanlocs, 'maplimits', 'minmax');
          topoplot(g.plottopo(index,:), g.chanlocs);
          if ~isempty(g.plottopotitle)
              title(g.plottopotitle{index}, 'interpreter', 'none');
          end
      end
      
      axes(tmpax);
  end
      
  for ind = 1:size(R,1)
      if ind == size(R,1) && strcmpi(g.plotmean, 'on') && size(R,1) > 1
           plot(times,R(ind,:), 'k', 'linewidth', 2);
      elseif ~isempty(g.colors),
           tmp = plot(times,R(ind,:), 'k'); 
           tmpcol = g.colors{mod(ind-1, length(g.colors))+1};
           if length(tmpcol) > 1, tmpstyle = tmpcol(2:end); tmpcol = tmpcol(1); else tmpstyle = '-'; end
           set(tmp, 'color', tmpcol, 'linestyle', tmpstyle); 
           
           if ~isempty(g.traceinfo)
               if ischar(g.traceinfo) && strcmpi(g.traceinfo, 'on')
                   set(tmp, 'ButtonDownFcn', [ 'disp(''Trace ' int2str(ind) ''');' ]);
               elseif iscell(g.traceinfo)
                   try
                       set(tmp, 'ButtonDownFcn', g.traceinfo{ind});
                   catch,
                       error('Trace info cell array does not contain the same number of element as trace in the graph')
                   end
               end
           end
           
           % change the line style when number of plots exceed number of colors in g.colors
           %lineStyles = {'-', '--',':','-.'};
           %set(tmp,'LineStyle',lineStyles{min(ceil(ind/length(g.colors)),length(lineStyles))});
          
           hold on;
      else plot(times,R(ind,:));
      end
  end
  
  % ordinate limits
  % ---------------
  if isempty(g.ylim), 
      yll = min(reshape(R, [1 prod(size(R))]));
      ylh = max(reshape(R, [1 prod(size(R))]));
      yll2 = yll - (ylh-yll)/10;
      ylh2 = ylh + (ylh-yll)/10;
      if ~isnan(yll), g.ylim = [yll2 ylh2]; end
  end
  if ~isempty(g.ylim) && length(g.ylim) == 2 
      if any(g.ylim)
          ylim(g.ylim);
      else
          ylim([0 1]);
          axis off;
          box off;
      end
  elseif ~isempty(g.ylim)
      yl = ylim;
      ylim([g.ylim yl(2)]);
  end
  yl = ylim; 

  % highlight regions
  % -----------------
  if ~isempty(g.maskarray)
      axsignif = highlight(ax, times, Rregions, g.highlightmode, g.xlabel);

      % replot data (above highlighted regions)
      % ---------
      axes(ax);
      for ind = 1:size(R,1)
          if ind == size(R,1) && strcmpi(g.plotmean, 'on') && size(R,1) > 1
               plot(times,R(ind,:), 'k', 'linewidth', 2);
          elseif ~isempty(g.colors),             
              tmp = plot(times,R(ind,:), 'k'); set(tmp, 'color', g.colors{mod(ind-1, length(g.colors))+1} ); hold on;
          else plot(times,R(ind,:));
          end
      end
      if strcmpi(g.highlightmode, 'bottom'), xlabel(''); set(ax, 'xtick', []); end
  end
  box on;
  
  ylim(yl);
  if strcmpi(g.logpval, 'on')
      set(gca, 'ytickmode', 'manual', 'yticklabel', round(10.^-get(gca, 'ytick')*1000)/1000, 'ydir', 'reverse');
  end
  
  % vertical lines
  % --------------
  hold on
  xl = xlim;
  if ~isnan(g.marktimes) % plot marked time
      for mt = g.marktimes(:)'
          plot([mt mt],[yl(1) yl(2)],'--k','LineWidth',g.linewidth);
      end
  end
  hold off
  if ~isempty(g.vert)
      for index = 1:length(g.vert)
          line([g.vert(index), g.vert(index)], [yl(1) yl(2)], 'linewidth', 1, 'color', 'm');
      end
  end
  xlim([times(1) times(end)]);

  % title and legend
  % ----------------
  if strcmpi(g.plotmode, 'topo') % plot in scalp array
      NAME_OFFSETX = 0.1;
      NAME_OFFSETY = 0.2;
      xx = xlim; xmin = xx(1); xdiff = xx(2)-xx(1); xpos = double(xmin+NAME_OFFSETX*xdiff);
      yy = ylim; ymax = yy(2); ydiff = yy(2)-yy(1); ypos = double(ymax-NAME_OFFSETY*ydiff);
      t=text(xpos, ypos,g.title);
      axis off;
      line([0 0], [yl(1) yl(2)], 'linewidth', 1, 'color', 'k');
      line([xl(1) xl(2)], [0 0], 'linewidth', 1, 'color', 'k');
      set(ax, 'userdata', { g.xlabel g.ylabel g.legend });
  else
      title(g.title, 'interpreter', 'none')
      if ~isempty(g.legend)
          hh = legend(g.legend(:), 'location', 'southeast');
          set(hh, 'unit', 'pixels', 'interpreter', 'none')
      end
      if isempty(g.maskarray)
          xlabel(g.xlabel);
      end
      ylabel(g.ylabel)
  end
  
% -----------------
% highlight regions
% -----------------
function axsignif = highlight(ax, times, regions, highlightmode, myxlabel);
color1 = [0.75 0.75 0.75];
color2 = [0 0 0];
yl  = ylim; 
yl(1) = yl(1)-max(abs(yl));
yl(2) = yl(2)+max(abs(yl));

if ~strcmpi(highlightmode, 'background')
    pos = get(ax, 'position');
    set(gca, 'xtick', []);
    axsignif = axes('position', [pos(1) pos(2)-pos(4)*0.05 pos(3) pos(4)*0.05 ]);
    plot(times, times, 'w');
    set(axsignif, 'ytick', []);
    yl2 = ylim;
    yl2(1) = yl2(1)-max(abs(yl2));
    yl2(2) = yl2(2)+max(abs(yl2));
    xlim([times(1) times(end)]);
    xlabel(myxlabel);
else
    axsignif = [];
    xlabel(myxlabel);
end

if ~isempty(regions)
    axes(ax);
    in_a_region = 0;
    for index=1:length(regions)
        if regions(index) && ~in_a_region
            tmpreg(1) = times(index);
            in_a_region = 1;
        end
        if (~regions(index) || index == length(regions)) && in_a_region
            tmpreg(2) = times(min(length(times), index));
            in_a_region = 0;
            if strcmpi(highlightmode, 'background')
                tmph = patch([tmpreg(1) tmpreg(2) tmpreg(2) tmpreg(1)], ...
                    [yl(1) yl(1) yl(2) yl(2)], color1); hold on;
                set(tmph, 'edgecolor', color1);
            else
                oldax = ax;
                axes(axsignif);
                tmph = patch([tmpreg(1) tmpreg(2) tmpreg(2) tmpreg(1)], ...
                    [yl2(1) yl2(1) yl2(2) yl2(2)], color2); hold on;
                set(tmph, 'edgecolor', color2);
                axes(oldax);
            end
        end
    end
    ylim(yl);
end
  

  function res = dims(array)
    res = min(ndims(array), max(size(array,2),size(array,3)));
