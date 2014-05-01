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

function std_chantopo(data, varargin)

pgroup = [];
pcond  = [];
pinter = [];
if nargin < 2
    help std_chantopo;
    return;
end;

opt = finputcheck( varargin, { 'ylim'        'real'    []              [];
                               'titles'      'cell'    []              cell(20,20);
                               'threshold'   'real'    []              NaN;
                               'chanlocs'    'struct'  []              struct('labels', {});
                               'groupstats'  'cell'    []              {};
                               'condstats'   'cell'    []              {};
                               'interstats'  'cell'    []              {};
                               'subplotpos'  'integer' []              [];
                               'topoplotopt' 'cell'    []              { 'style', 'both' };
                               'binarypval'  'string'  { 'on','off' }  'on';
                               'datatype'    'string'  { 'ersp','itc','erp','spec' }    'erp';
                               'caxis'       'real'    []              [] }, 'std_chantopo', 'ignore'); %, 'ignore');
if isstr(opt), error(opt); end;
if ~isempty(opt.ylim), opt.caxis = opt.ylim; end;
if isnan(opt.threshold), opt.binarypval = 'off'; end;
if strcmpi(opt.binarypval, 'on'), opt.ptopoopt = { 'style' 'blank' }; else opt.ptopoopt = opt.topoplotopt; end;

% remove empty entries
datapresent = ~cellfun(@isempty, data);
for c = size(data,1):-1:1, if sum(datapresent(c,:)) == 0, data(c,:) = []; opt.titles(c,:) = []; if ~isempty(opt.groupstats), opt.groupstats(c) = []; end; end; end;
for g = size(data,2):-1:1, if sum(datapresent(:,g)) == 0, data(:,g) = []; opt.titles(:,g) = []; if ~isempty(opt.condstats ), opt.condstats( g) = []; end; end; end;

nc = size(data,1);
ng = size(data,2);
if nc >= ng, opt.transpose = 'on';
else         opt.transpose = 'off';
end;

% plotting paramters
% ------------------
if ng > 1 & ~isempty(opt.groupstats), addc = 1; else addc = 0; end;
if nc > 1 & ~isempty(opt.condstats ), addr = 1; else addr = 0; end;
if ~isempty(opt.subplotpos), 
     if strcmpi(opt.transpose, 'on'), opt.subplotpos = opt.subplotpos([2 1 4 3]); end;
     addr = opt.subplotpos(1);
     addc = opt.subplotpos(2);
     posr = opt.subplotpos(4);
     posc = opt.subplotpos(3);
else posr = 0;
     posc = 0;
end;

% compute significance mask
% -------------------------
if ~isempty(opt.interstats), pinter = opt.interstats{3}; end;

if ~isnan(opt.threshold(1)) && ( ~isempty(opt.groupstats) || ~isempty(opt.condstats) )    
    pcondplot  = opt.condstats;
    pgroupplot = opt.groupstats;
    pinterplot = pinter;
    maxplot = 1;
else
    for ind = 1:length(opt.condstats),  pcondplot{ind}  = -log10(opt.condstats{ind}); end;
    for ind = 1:length(opt.groupstats), pgroupplot{ind} = -log10(opt.groupstats{ind}); end;
    if ~isempty(pinter), pinterplot = -log10(pinter); end;
    maxplot = 3;
end;

% adjust figure size
% ------------------
if isempty(opt.subplotpos)
    fig = figure('color', 'w');
    pos = get(fig, 'position');
    set(fig, 'position', [ pos(1)+15 pos(2)+15 pos(3)/2.5*(nc+addr), pos(4)/2*(ng+addc) ]);
    pos = get(fig, 'position');
    if strcmpi(opt.transpose, 'off'), set(gcf, 'position', [ pos(1) pos(2) pos(4) pos(3)]);
    else                              set(gcf, 'position', pos);
    end;
end

% topoplot
% --------
tmpc = [inf -inf];
for c = 1:nc
    for g = 1:ng
        hdl(c,g) = mysubplot(nc+addr, ng+addc, g + posr + (c-1+posc)*(ng+addc), opt.transpose);
        if ~isempty(data{c,g})
            tmpplot = double(mean(data{c,g},3));
            if ~all(isnan(tmpplot))
                if ~isreal(tmpplot), error('This function cannot plot complex values'); end;
                topoplot( tmpplot, opt.chanlocs, opt.topoplotopt{:});
                if isempty(opt.caxis)
                    tmpc = [ min(min(tmpplot), tmpc(1)) max(max(tmpplot), tmpc(2)) ];
                else 
                    caxis(opt.caxis);
                end;
                title(opt.titles{c,g}, 'interpreter', 'none'); 
            else
                axis off;
            end;
        else
            axis off;
        end;

        % statistics accross groups
        % -------------------------
        if g == ng & ng > 1 & ~isempty(opt.groupstats)
            hdl(c,g+1) = mysubplot(nc+addr, ng+addc, g + posr + 1 + (c-1+posc)*(ng+addc), opt.transpose);
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
    if ~isempty(opt.condstats) && nc > 1
        hdl(nc+1,g) = mysubplot(nc+addr, ng+addc, g + posr + (c+posc)*(ng+addc), opt.transpose);
        topoplot( pcondplot{g}, opt.chanlocs, opt.ptopoopt{:});
        title(opt.titles{nc+1,g}); 
        caxis([-maxplot maxplot]);
    end;
end;

% statistics accross group and conditions
% ---------------------------------------
if ~isempty(opt.condstats) && ~isempty(opt.groupstats) && ng > 1 && nc > 1
    hdl(nc+1,ng+1) = mysubplot(nc+addr, ng+addc, g + posr + 1 + (c+posc)*(ng+addc), opt.transpose);
    topoplot( pinterplot, opt.chanlocs, opt.ptopoopt{:});
    title(opt.titles{nc+1,ng+1}); 
    caxis([-maxplot maxplot]);
end;    

% color bars
% ----------
axes(hdl(nc,ng)); 
cbar_standard(opt.datatype, ng);
if isnan(opt.threshold(1)) && (nc ~= size(hdl,1) || ng ~= size(hdl,2))
    axes(hdl(end,end));
    cbar_signif(ng, maxplot);
end;

% remove axis labels
% ------------------
for c = 1:size(hdl,1)
    for g = 1:size(hdl,2)
        if g ~= 1 && size(hdl,2) ~=1, ylabel(''); end;
        if c ~= size(hdl,1) && size(hdl,1) ~= 1, xlabel(''); end;
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

