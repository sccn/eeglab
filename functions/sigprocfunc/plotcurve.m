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
%
% Authors: Arnaud Delorme, 2004, Bhaktivedanta Institute

% Copyright (C) 2004  Arnaud Delorme
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

function plotcurve( times, R, varargin);

	if nargin < 2
        help plotcurve;
        return;
	end;
	g = finputcheck( varargin, { 'maskarray'     'real'    []       [];
                             'val2mask'      'real'    []           R;
                             'highlightmode' 'string'  { 'background' 'bottom' } 'background';
                             'plotmean'      'string'  { 'on' 'off' }            'on';
                             'title'         'string'  []                        '';
                             'xlabel'        'string'  []                        '';
                             'ylabel'        'string'  []                        '';
                             'legend'        'cell'    []                        {};
                             'ylim'          'real'    []                        [];
                             'vert'          'real'    []                        [];
                             'linewidth'     'real'    []                        2;
                             'marktimes'     'real'    []                        [] });
   if isstr(g), error(g); end;

   % regions of significance
   % -----------------------
   if ~isempty(g.maskarray)
      Rregions = ones(size(g.val2mask));
      switch dims(g.maskarray)
           case 3, Rregions  (find(g.val2mask > g.maskarray(:,:,1) & (g.val2mask < g.maskarray(:,:,2)))) = 0;
           case 2, if size(g.val2mask,2) == size(g.maskarray,2)
                       Rregions  (find(g.val2mask < g.maskarray)) = 0;
                   else
                       Rregions(find((g.val2mask > repmat(g.maskarray(:,1),[1 length(times)])) ...
                         & (g.val2mask < repmat(g.maskarray(:,2),[1 length(times)])))) = 0;
                   end;
           case 1, Rregions  (find(g.val2mask < repmat(g.maskarray(:),[1 size(g.val2mask,2)]))) = 0;
       end; 
       Rregions = sum(Rregions,1);
   else 
       Rregions = [];
   end

  % plotting
  % --------
  if strcmpi(g.plotmean, 'on')
      R       = mean(R,1);
  end;
  plot(times,R);
  if ~isempty(g.ylim) && length(g.ylim) == 2
      ylim(g.ylim);
  elseif ~isempty(g.ylim)
      yl = ylim;
      ylim([g.ylim yl(2)]);
  end
  if ~isempty(g.maskarray)
      highlight(gca, times, Rregions, g.highlightmode, g.xlabel);
  else
      xlabel(g.xlabel);
  end;
  if strcmpi(g.highlightmode, 'background')
      plot(times,R(:,:)); % replot
  end;

  % vertical lines
  % --------------
  hold on
  yl = ylim;
  if times(1) < 0
      plot([0 0],[yl(1) yl(2)],'--m','LineWidth',g.linewidth);
  end;
  if ~isnan(g.marktimes) % plot marked time
      for mt = g.marktimes(:)'
          plot([mt mt],[yl(1) yl(2)],'--k','LineWidth',g.linewidth);
      end
  end
  hold off
  if ~isempty(g.vert)
      for index = 1:length(g.vert)
          line([g.vert(index), g.vert(index)], [yl(1) yl(2)], 'linewidth', 1, 'color', 'm');
      end;
  end;
  xlim([times(1) times(end)]);

  % title and legend
  % ----------------
  title(g.title)
  ylabel(g.ylabel)
  if ~isempty(g.legend)
      legend(g.legend);
  end;
  
function highlight(ax, times, regions, highlightmode, myxlabel);
color1 = [0.75 0.75 0.75];
color2 = [0 0 0];
yl  = ylim;

if ~strcmpi(highlightmode, 'background')
    pos = get(gca, 'position');
    set(gca, 'xtick', []);
    axsignif = axes('position', [pos(1) pos(2)-pos(4)*0.05 pos(3) pos(4)*0.05 ]);
    plot(times, times, 'w');
    set(axsignif, 'ytick', []);
    yl2 = ylim;
    xlim([times(1) times(end)]);
    xlabel(myxlabel);
else
    xlabel(myxlabel);
end;

if ~isempty(regions)
    axes(ax);
    in_a_region = 0;
    for index=1:length(regions)
        if regions(index) & ~in_a_region
            tmpreg(1) = times(index);
            in_a_region = 1;
        end;
        if ~regions(index) & in_a_region
            tmpreg(2) = times(index);
            in_a_region = 0;
            if strcmpi(highlightmode, 'background')
                tmph = patch([tmpreg(1) tmpreg(2) tmpreg(2) tmpreg(1)], ...
                    [yl(1) yl(1) yl(2) yl(2)], color1); hold on;
                set(tmph, 'edgecolor', color1);
            else
                oldax = gca;
                axes(axsignif);
                tmph = patch([tmpreg(1) tmpreg(2) tmpreg(2) tmpreg(1)], ...
                    [yl2(1) yl2(1) yl2(2) yl2(2)], color2); hold on;
                set(tmph, 'edgecolor', color2);
                axes(oldax);
            end;
        end;
    end;
end;
  
  function res = dims(array)
    res = min(ndims(array), max(size(array,2),size(array,3)));
