% std_chantopo() - plot ERP/spectral/ERSP topoplot at a specific
%                  latency/frequency. 
% Usage:
%          >> std_chantopo( data, 'key', 'val', ...)
% Inputs:
%  data  -  [cell array] mean data for each subject group and/or data
%           condition. These arrays are usually returned by function
%           std_erspplot and std_erpplot. For example
%
%           >> data = { [1x64x12] [1x64x12 }; % 2 groups of 12 subjects, 64 channels
%           >> std_chantopo(data, 'chanlocs', 'chanlocfile.txt');
%
% Scalp map plotting option (mandatory):
%  'chanlocs'    - [struct] channel location structure
%
% Other scalp map plotting options:
%  'chanlocs'    - [struct] channel location structure
%  'topoplotopt' - [cell] topoplot options. Default is { 'style', 'both', 
%                  'shading', 'interp' }. See topoplot help for details.
%  'ylim'        - [min max] ordinate limits for ERP and spectrum plots
%                  {default: all available data}
%  'caxis'       - [min max] same as above
%
% Optional display parameters:
%  'datatype'    - ['erp'|'spec'] data type {default: 'erp'}
%  'titles'      - [cell array of string] titles for each of the subplots. 
%                  { default: none}
%  'subplotpos'  - [addr addc posr posc] perform ploting in existing figure.
%                  Add "addr" rows, "addc" columns and plot the scalp
%                  topographies starting at position (posr,posc).
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
% Author: Arnaud Delorme, CERCO, CNRS, 2006-
%
% See also: pop_erspparams(), pop_erpparams(), pop_specparams(), statcond()

% Copyright (C) 2006 Arnaud Delorme
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

function std_chantopo(data, varargin)

pgroup = [];
pcond  = [];
pinter = [];
if nargin < 2
    help std_chantopo;
    return;
end

opt = finputcheck( varargin, { 'ylim'        'real'    []              [];
                               'titles'      'cell'    []              cell(20,20);
                               'threshold'   'real'    []              NaN;
                               'chanlocs'    'struct'  []              struct('labels', {});
                               'groupstats'  'cell'    []              {};
                               'condstats'   'cell'    []              {};
                               'interstats'  'cell'    []              {};
                               'effect'      'string' { 'main','marginal' }   'marginal';
                               'subplotpos'  'integer' []              [];
                               'topoplotopt' 'cell'    []              { 'style', 'both' };
                               'binarypval'  'string'  { 'on','off' }  'on';
                               'datatype'    'string'  { 'ersp','itc','erp','spec' }    'erp';
                               'caxis'       'real'    []              [] }, 'std_chantopo', 'ignore'); %, 'ignore');
if ischar(opt), error(opt); end
if ~isempty(opt.ylim), opt.caxis = opt.ylim; end
if isnan(opt.threshold), opt.binarypval = 'off'; end
if strcmpi(opt.binarypval, 'on'), opt.ptopoopt = { 'style' 'blank' }; else opt.ptopoopt = opt.topoplotopt; end

% remove empty entries
datapresent = ~cellfun(@isempty, data);
for c = size(data,1):-1:1, if sum(datapresent(c,:)) == 0, data(c,:) = []; opt.titles(c,:) = []; if ~isempty(opt.groupstats), opt.groupstats(c) = []; end; end; end
for g = size(data,2):-1:1, if sum(datapresent(:,g)) == 0, data(:,g) = []; opt.titles(:,g) = []; if ~isempty(opt.condstats ), opt.condstats( g) = []; end; end; end

nc = size(data,1);
ng = size(data,2);
if nc >= ng, opt.transpose = 'on';
else         opt.transpose = 'off';
end

% plotting paramters
% ------------------
if ng > 1 && ~isempty(opt.groupstats), addc = 1; else addc = 0; end
if nc > 1 && ~isempty(opt.condstats ), addr = 1; else addr = 0; end
if ~isempty(opt.subplotpos), 
     if strcmpi(opt.transpose, 'on'), opt.subplotpos = opt.subplotpos([2 1 4 3]); end
     addr = opt.subplotpos(1);
     addc = opt.subplotpos(2);
     posr = opt.subplotpos(4);
     posc = opt.subplotpos(3);
else posr = 0;
     posc = 0;
end

% compute significance mask
% -------------------------
pinterplot = {};
if strcmpi(opt.effect, 'marginal') || ng == 1 || nc == 1
    if ~isnan(opt.threshold) && ( ~isempty(opt.groupstats) || ~isempty(opt.condstats) )    
        pcondplot  = opt.condstats;
        pgroupplot = opt.groupstats;
        maxplot = 1;
    else
        for ind = 1:length(opt.condstats),  pcondplot{ind}  = -log10(opt.condstats{ind}); end
        for ind = 1:length(opt.groupstats), pgroupplot{ind} = -log10(opt.groupstats{ind}); end
        maxplot = 3;
    end
elseif strcmpi(opt.effect, 'main') && ~isempty(opt.interstats)
    if ~isnan(opt.threshold) && ( ~isempty(opt.groupstats) || ~isempty(opt.condstats) )    
        pcondplot  = { opt.interstats{1} };
        pgroupplot = { opt.interstats{2} };
        pinterplot = opt.interstats{3};
        maxplot = 1;
    else
        if ~isempty(opt.interstats{1}), pcondplot  = { -log10(opt.interstats{1}) }; end
        if ~isempty(opt.interstats{2}), pgroupplot = { -log10(opt.interstats{2}) }; end
        if ~isempty(opt.interstats{3}), pinterplot = -log10(opt.interstats{3}); end
        maxplot = 3;
    end
end

% adjust figure size
% ------------------
if isempty(opt.subplotpos)
    fig = figure('color', 'w');
    pos = get(fig, 'position');
    set(fig, 'position', [ pos(1)+15 pos(2)+15 pos(3)/2.5*(nc+addr), pos(4)/2*(ng+addc) ]);
    pos = get(fig, 'position');
    if strcmpi(opt.transpose, 'off'), set(gcf, 'position', [ pos(1) pos(2) pos(4) pos(3)]);
    else                              set(gcf, 'position', pos);
    end
end

% topoplot
% --------
tmpc = [inf -inf];
for c = 1:nc
    for g = 1:ng
        hdl(c,g) = mysubplot(nc+addr, ng+addc, c, g, opt.transpose);
        if ~isempty(data{c,g})
            tmpplot = double(mean(data{c,g},3));
            if ~isreal(tmpplot(1)), tmpplot = abs(tmpplot); end % comes second for processing single trials
            if ~all(isnan(tmpplot))
                if ~isreal(tmpplot), error('This function cannot plot complex values'); end
                topoplot( tmpplot, opt.chanlocs, opt.topoplotopt{:});
                if isempty(opt.caxis)
                    tmpc = [ min(min(tmpplot), tmpc(1)) max(max(tmpplot), tmpc(2)) ];
                else 
                    caxis(opt.caxis);
                end
                title(opt.titles{c,g}, 'interpreter', 'none'); 
            else
                axis off;
            end
        else
            axis off;
        end

        % statistics accross groups
        % -------------------------
        if strcmpi(opt.effect, 'marginal') || (strcmpi(opt.effect, 'main') && c == 1)
            if g == ng && ng > 1 && ~isempty(opt.groupstats)
                if strcmpi(opt.effect, 'main') && nc>1, centerc = nc/2-0.5; else centerc = 0; end
                hdl(c,g+1) = mysubplot(nc+addr, ng+addc, c+centerc, ng + 1, opt.transpose);
                pgroupplot{c}(pgroupplot{c} <0) = 0;
                topoplot( pgroupplot{c}, opt.chanlocs, opt.ptopoopt{:});
                title(opt.titles{c,g+1});
                caxis([-maxplot maxplot]);
            end
        end
    end
end
        
% color scale
% -----------
if isempty(opt.caxis)
    for c = 1:nc
        for g = 1:ng
            axes(hdl(c,g));
            caxis(tmpc);
        end
    end
end

for g = 1:ng
    % statistics accross conditions
    % -----------------------------
    if strcmpi(opt.effect, 'marginal') || (strcmpi(opt.effect, 'main') && g == 1)
        if ~isempty(opt.condstats) && nc > 1
            if strcmpi(opt.effect, 'main') && ng>1, centerg = ng/2-0.5; else centerg = 0; end
            hdl(nc+1,g) = mysubplot(nc+addr, ng+addc, nc+addr, g+centerg, opt.transpose);
            pcondplot{g}(pcondplot{g} < 0) = 0;
            topoplot( pcondplot{g}, opt.chanlocs, opt.ptopoopt{:});
            title(opt.titles{nc+1,g});
            caxis([-maxplot maxplot]);
        end
    end
end

% statistics accross group and conditions
% ---------------------------------------
if ~isempty(opt.condstats) && ~isempty(opt.groupstats) && ng > 1 && nc > 1 && ~isempty(pinterplot)
    hdl(nc+1,ng+1) = mysubplot(nc+addr, ng+addc, nc+addr, ng+1, opt.transpose);
    pinterplot(pinterplot<0) = 0;
    topoplot( pinterplot, opt.chanlocs, opt.ptopoopt{:});
    title(opt.titles{nc+1,ng+1}); 
    caxis([-maxplot maxplot]);
end    

% color bars
% ----------
if isnan(opt.threshold(1)) && (nc ~= size(hdl,1) || ng ~= size(hdl,2))
    try
        axes(hdl(end,end));
        cbar_signif(ng, maxplot);
    catch
        cbar_signif(ng, maxplot);
    end
end
axes(hdl(nc,ng)); 
cbar_standard(opt.datatype, ng);

% remove axis labels
% ------------------
for c = 1:size(hdl,1)
    for g = 1:size(hdl,2)
        if g ~= 1 && size(hdl,2) ~=1, ylabel(''); end
        if c ~= size(hdl,1) && size(hdl,1) ~= 1, xlabel(''); end
    end
end

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
    end

% colorbar for significance
% -------------------------
function cbar_signif(ng, maxplot);
    % Retrieving Defaults
    icadefs;
    
    pos = get(gca, 'position');
    tmpc = caxis;
    fact = fastif(ng == 1, 40, 20);
    tmp = axes('position', [ pos(1)+pos(3)+pos(3)/fact pos(2) pos(3)/fact pos(4) ]);  
    map = colormap(DEFAULT_COLORMAP);
    n = size(map,1);
    cols = [ceil(n/2):n]';
    image([0 1],linspace(0,maxplot,length(cols)),[cols cols]);
    %cbar(tmp, 0, tmpc, 5);
    tick = linspace(0, maxplot, maxplot+1);
    set(gca, 'ytickmode', 'manual', 'YAxisLocation', 'right', 'xtick', [], ...
        'ytick', tick, 'yticklabel', round(10.^-tick*1000)/1000);
    xlabel('');
    colormap(DEFAULT_COLORMAP);

% mysubplot (allow to transpose if necessary)
% -------------------------------------------
function hdl = mysubplot(nr,nc,r,c,subplottype)

    cmargin = 0.2/nc;
    rmargin = 0.2/nr;
    if strcmpi(subplottype, 'transpose') || strcmpi(subplottype, 'on'),   hdl = subplot('position',[(r-1)/nr+rmargin (nc-c)/nc+cmargin 1/nr-2*rmargin 1/nc-2*cmargin]);
    elseif strcmpi(subplottype, 'normal') || strcmpi(subplottype, 'off'), hdl = subplot('position',[(c-1)/nc+cmargin (nr-r)/nr+rmargin 1/nc-2*cmargin 1/nr-2*rmargin]);
    elseif strcmpi(subplottype, 'noplot'), hdl = gca;
    else error('Unknown subplot type');
    end
