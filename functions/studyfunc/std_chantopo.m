% std_chantopo() - plot ERP/spectral/ERSP topoplot at a specific
%                  latency/frequency. 
% Usage:
%          >> std_chantopo( data, 'key', 'val', ...)
% Inputs:
%  data  -  [cell array] mean data for each subject group and/or data
%           condition. For example, to plot mean ERPs from a STUDY 
%           for epochs of 800 frames in two conditions from three groups 
%           of 12 subjects:
%
%           >> data = { [800x12] [800x12] [800x12];... % 3 groups, cond 1
%                       [800x12] [800x12] [800x12] };  % 3 groups, cond 2
%           >> std_chantopo(erp_ms,data);
%
%           By default, parametric statistics are computed across subjects 
%           in the three groups. (group,condition) ERP averages are plotted. 
%           See below and >> help statcond 
%           for more information about the statistical computations.
%
% Optional display parameters:
%  'datatype'    - ['erp'|'spec'] data type {default: 'erp'}
%  'titles'      - [cell array of string] titles for each of the subplots. 
%                  { default: none}
%
% Statistics options:
%  'groupstats'  - [cell] One p-value array per group {default: {}}
%  'condstats'   - [cell] One p-value array per condition {default: {}}
%  'interstats'  - [cell] Interaction p-value arrays {default: {}}
%  'threshold'   - [NaN|real<<1] Significance threshold. NaN -> plot the 
%                  p-values themselves on a different figure. When possible, 
%                  significance regions are indicated below the data.
%                  {default: NaN}
%  'binarypval'  - ['on'|'off'] if a threshold is set, show only significant
%                  channels as red dots. Default is 'off'.
%
% Curve plotting options:
% 'ylim'         - [min max] ordinate limits for ERP and spectrum plots
%                  {default: all available data}
% 'caxis'        - [min max] same as above
%
% Scalp map plotting options:
%  'chanlocs'    - [struct] channel location structure
%
% Author: Arnaud Delorme, CERCO, CNRS, 2006-
%
% See also: pop_erspparams(), pop_erpparams(), pop_specparams(), statcond()

% $Log: not supported by cvs2svn $
% Revision 1.12  2010/02/09 06:07:27  arno
% Fixed new title problem and implemented 3-level significance
%
% Revision 1.11  2010/02/06 05:47:52  arno
% New titles for figures
%
% Revision 1.10  2010/02/01 19:00:48  arno
% New binary pval option
%
% Revision 1.9  2009/11/23 22:13:20  arno
% Fixed the caxis problem
%
% Revision 1.8  2009/10/10 01:31:23  arno
% fix scalp maps plotting for ERP (single subjects)
%
% Revision 1.7  2009/08/29 04:24:56  arno
% new statistics
%
% Revision 1.5  2009/08/29 00:38:31  arno
% move all statistics to std_stat
%
% Revision 1.4  2009/08/11 00:22:59  arno
% fix bootstrap problem
%
% Revision 1.3  2009/05/31 02:22:10  arno
% Adding FDR and bootstrap to all STUDY functions
%
% Revision 1.2  2008/03/30 12:04:11  arno
% fix minor problem (title...)
%
% Revision 1.1  2007/01/26 18:08:17  arno
% Initial revision
%
% Revision 1.24  2006/11/23 00:26:22  arno
% cosmetic change
%
% Revision 1.23  2006/11/22 20:22:38  arno
% filter if necessary
%
% Revision 1.22  2006/11/22 19:34:22  arno
% empty plot
%
% Revision 1.21  2006/11/16 00:07:31  arno
% fix title for topoplot
%
% Revision 1.20  2006/11/15 21:58:49  arno
% plotting titles
%
% Revision 1.19  2006/11/09 23:55:18  arno
% fix last change
%
% Revision 1.18  2006/11/09 23:27:11  arno
% figure titles
%
% Revision 1.17  2006/11/09 23:21:21  arno
% add component name
%
% Revision 1.16  2006/11/09 22:04:39  arno
% ERSP plotting
%
% Revision 1.15  2006/11/03 03:01:17  arno
% allow plotting specific time-freq point
%
% Revision 1.14  2006/11/02 22:14:15  arno
% g. -> opt.
%
% Revision 1.13  2006/11/02 22:13:04  arno
% Limit to ERSP
%
% Revision 1.12  2006/11/02 21:59:18  arno
% input option
%
% Revision 1.11  2006/11/02 21:54:51  arno
% same
%
% Revision 1.10  2006/11/02 21:53:06  arno
% condstats -> condstat
%
% Revision 1.9  2006/10/10 23:50:54  scott
% replaced ?? defaults with defaults from finputcheck()
%
% Revision 1.8  2006/10/09 23:51:45  scott
% some more help edits
%
% Revision 1.7  2006/10/04 23:55:34  toby
% Bug fix courtesy Bas de Kruif
%
% Revision 1.6  2006/10/03 21:46:11  scott
% edit help msg -- ?? remain... -sm
%
% Revision 1.5  2006/10/03 16:29:27  scott
% edit
%
% Revision 1.4  2006/10/03 16:24:12  scott
% help message eidts.  ARNO - SEE ??s   -sm
%
% Revision 1.3  2006/10/02 11:41:13  arno
% wrote documentation
%
% Copyright (C) 2006 Arnaud Delorme
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

function std_chantopo(data, varargin)

pgroup = [];
pcond  = [];
pinter = [];
if nargin < 2
    help std_chantopo;
    return;
end;

opt = finputcheck( varargin, { 'ylim'        'real'   []              [];
                               'titles'      'cell'   []              cell(20,20);
                               'threshold'   'real'   []              NaN;
                               'chanlocs'    'struct' []              struct('labels', {});
                               'groupstats'  'cell'   []              {};
                               'condstats'   'cell'   []              {};
                               'interstats'  'cell'   []              {};
                               'topoplotopt' 'cell'   []              { 'style', 'map', 'shading', 'interp' };
                               'binarypval'  'string' { 'on' 'off' }  'off';
                               'datatype'    'string' { 'ersp' 'itc' 'erp' 'spec' }    'erp';
                               'caxis'       'real'   []              [] }, 'std_chantopo', 'ignore'); %, 'ignore');
if isstr(opt), error(opt); end;
if ~isempty(opt.ylim), opt.caxis = opt.ylim; end;
if strcmpi(opt.binarypval, 'on'), opt.ptopoopt = { 'style' 'blank' }; else opt.ptopoopt = opt.topoplotopt; end;
if isempty(opt.titles), opt.titles = cell(10,10); opt.titles(:) = { '' }; end;

nc = size(data,1);
ng = size(data,2);
if nc >= ng, opt.transpose = 'on';
else         opt.transpose = 'off';
end;

% plotting paramters
% ------------------
if ng > 1 & ~isempty(opt.groupstats), addc = 1; else addc = 0; end;
if nc > 1 & ~isempty(opt.condstats ), addr = 1; else addr = 0; end;

% compute significance mask
% -------------------------
if ~isempty(opt.interstats), pinter = opt.interstats{3}; end;

if ~isnan(opt.threshold) && ( ~isempty(opt.groupstats) || ~isempty(opt.condstats) )    
    pcondplot  = opt.condstats;
    pgroupplot = opt.groupstats;
    maxplot = 1;
else
    warning off;
    for ind = 1:length(opt.condstats),  pcondplot{ind}  = -log10(opt.condstats{ind}); end;
    for ind = 1:length(opt.groupstats), pgroupplot{ind} = -log10(opt.groupstats{ind}); end;
    if ~isempty(pinter), pinterplot = -log10(pinter); end;
    maxplot = 3;
    warning on;
end;

% adjust figure size
% ------------------
figure('color', 'w');
pos = get(gcf, 'position');
basewinsize = 200/max(nc,ng)*3;
pos(3) = 200*(ng+addc);
pos(4) = 200*(nc+addr);
if strcmpi(opt.transpose, 'on'), set(gcf, 'position', [ pos(1) pos(2) pos(4) pos(3)]);
else                             set(gcf, 'position', pos);
end;

% topoplot
% --------
tmpc = [inf -inf];
for c = 1:nc
    for g = 1:ng
        hdl(c,g) = mysubplot(nc+addr, ng+addc, g + (c-1)*(ng+addc), opt.transpose);
        if ~isempty(data{c,g})
            tmpplot = double(mean(data{c,g},3));
            topoplot( tmpplot, opt.chanlocs, opt.topoplotopt{:});
            if isempty(opt.caxis)
                tmpc = [ min(min(tmpplot), tmpc(1)) max(max(tmpplot), tmpc(2)) ];
            else 
                caxis(opt.caxis);
            end;
            title(opt.titles{c,g}); 
        else
            axis off;
        end;

        % statistics accross groups
        % -------------------------
        if g == ng & ng > 1 & ~isempty(opt.groupstats)
            hdl(c,g+1) = mysubplot(nc+addr, ng+addc, g + 1 + (c-1)*(ng+addc), opt.transpose);
            topoplot( pgroupplot{c}, opt.chanlocs, opt.ptopoopt{:});
            title(opt.titles{c,g+1}); 
            caxis([-maxplot maxplot]);
        end;
    end;
end;
        
% color scale
% -----------
if isempty(opt.caxis)
    for c = 1:nc
        for g = 1:ng
            axes(hdl(c,g));
            caxis(tmpc);
        end;
    end;
end;

for g = 1:ng
    % statistics accross conditions
    % -----------------------------
    if ~isempty(opt.condstats) & nc > 1
        hdl(nc+1,g) = mysubplot(nc+addr, ng+addc, g + c*(ng+addc), opt.transpose);
        topoplot( pcondplot{g}, opt.chanlocs, opt.ptopoopt{:});
        title(opt.titles{nc+1,g}); 
        caxis([-maxplot maxplot]);
    end;
end;

% statistics accross group and conditions
% ---------------------------------------
if ~isempty(opt.condstats) & ~isempty(opt.groupstats) & ng > 1 & nc > 1
    hdl(nc+1,ng+1) = mysubplot(nc+addr, ng+addc, g + 1 + c*(ng+addr), opt.transpose);
    topoplot( pinterplot, opt.chanlocs, opt.ptopoopt{:});
    title(opt.titles{nc+1,ng+1}); 
    caxis([-maxplot maxplot]);
end;    

% color bars
% ----------
axes(hdl(nc,ng)); 
cbar_standard(opt.datatype, ng);
if nc ~= size(hdl,1) | ng ~= size(hdl,2)
    axes(hdl(end,end));
    cbar_signif(ng, maxplot);
end;

% remove axis labels
% ------------------
for c = 1:size(hdl,1)
    for g = 1:size(hdl,2)
        if g ~= 1 & size(hdl,2) ~=1, ylabel(''); end;
        if c ~= size(hdl,1) & size(hdl,1) ~= 1, xlabel(''); end;
    end;
end;

% colorbar for ERSP and scalp plot
% --------------------------------
function cbar_standard(datatype, ng);
    pos = get(gca, 'position');
    tmpc = caxis;
    fact = fastif(ng == 1, 40, 20);
    tmp = axes('position', [ pos(1)+pos(3)+pos(3)/fact pos(2) pos(3)/fact pos(4) ]);  
    set(gca, 'unit', 'normalized');
    if strcmpi(datatype, 'itc')
         cbar(tmp, 0, tmpc, 10); ylim([0.5 1]);
    else cbar(tmp, 0, tmpc, 5);
    end;

% colorbar for significance
% -------------------------
function cbar_signif(ng, maxplot);
    pos = get(gca, 'position');
    tmpc = caxis;
    fact = fastif(ng == 1, 40, 20);
    tmp = axes('position', [ pos(1)+pos(3)+pos(3)/fact pos(2) pos(3)/fact pos(4) ]);  
    map = colormap;
    n = size(map,1);
    cols = [ceil(n/2):n]';
    image([0 1],linspace(0,maxplot,length(cols)),[cols cols]);
    %cbar(tmp, 0, tmpc, 5);
    tick = linspace(0, maxplot, maxplot+1);
    set(gca, 'ytickmode', 'manual', 'YAxisLocation', 'right', 'xtick', [], ...
        'ytick', tick, 'yticklabel', round(10.^-tick*1000)/1000);
    xlabel('');

% mysubplot (allow to transpose if necessary)
% -------------------------------------------
function hdl = mysubplot(nr,nc,ind,transp);

    r = ceil(ind/nc);
    c = ind -(r-1)*nc;
    if strcmpi(transp, 'on'), hdl = subplot(nc,nr,(c-1)*nr+r);
    else                      hdl = subplot(nr,nc,(r-1)*nc+c);
    end;

