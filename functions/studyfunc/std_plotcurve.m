% std_plotcurve() - plot ERP or spectral traces for a STUDY component 
%                   or channel cluster 
% Usage:
%          >> std_plotcurve( axvals, data, 'key', 'val', ...)
% Inputs:
%  axvals - [vector or cell array] axis values for the data. 
%  data  -  [cell array] mean data for each subject group and/or data
%           condition. For example, to plot mean ERPs from a STUDY 
%           for epochs of 800 frames in two conditions from three groups 
%           of 12 subjects:
%
%           >> data = { [800x12] [800x12] [800x12];... % 3 groups, cond 1
%                       [800x12] [800x12] [800x12] };  % 3 groups, cond 2
%           >> std_plotcurve(erp_ms,data);
%
%           By default, parametric statistics are computed across subjects 
%           in the three groups. (group,condition) ERP averages are plotted. 
%           See below and >> help statcond 
%           for more information about the statistical computations. For 
%           plotting multiple channels, use the second dimension. For
%           example data = { [800x64x12] [800x64x12] } for 12 subjects, 
%           64 channels and 800 data points. The 'chanlocs' option must be
%           used as well to specify channel positions.
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
%
% Curve plotting options (ERP and spectrum):
%  'plotgroups'  - ['together'|'apart'] 'together' -> plot mean results 
%                  for subject groups in the same figure panel in different 
%                  colors. 'apart' -> plot group results on different figure
%                  panels {default: 'apart'}
%  'plotconditions' - ['together'|'apart'] 'together' -> plot mean results 
%                  for data conditions on the same figure panel in
%                  different 
%                  colors. 'apart' -> plot conditions on different figure
%                  panel. Note: 'plotgroups' and 'plotconditions' arguments 
%                  cannot both be 'together' {default: 'apart'}
%  'legend'      - ['on'|'off'] turn plot legend on/off {default: 'off'}
%  'plotdiff'    - ['on'|'off'] plot difference between two groups
%                  or conditions plotted together. 
%  'plotstderr'  - ['on'|'off'|'diff'|'nocurve'|'diffnocurve'] plots in 
%                  a surface indicating the standard error. 'diff' only 
%                  does it for the difference (requires 'plotdiff' 'on' 
%                  above). 'nocurve' does not plot the mean. This functionality
%                  does not work for all data configuration {default: 'off'}
%  'figure'      - ['on'|'off'] creates a new figure ('on'). The 'off' mode
%                  plots all of the groups and conditions on the same pannel.
% 'plotsubjects' - ['on'|'off'] overplot traces for individual components
%                  or channels {default: 'off'}
% 'singlesubject' - ['on'|'off'] set to 'on' to plot single subject.
%                  {default: 'off'}
% 'ylim'         - [min max] ordinate limits for ERP and spectrum plots
%                  {default: all available data}
%
% Scalp map plotting options:
%  'chanlocs'    - [struct] channel locations structure
%
% Author: Arnaud Delorme, CERCO, CNRS, 2006-
%
% See also: pop_erspparams(), pop_erpparams(), pop_specparams(), statcond()

function std_plotcurve(allx, data, varargin) 

pgroup = [];
pcond  = [];
pinter = [];
if nargin < 2
    help std_plotcurve;
    return;
end;

opt = finputcheck( varargin, { 'ylim'          'real'   []              [];
                               'filter'        'real'   []              [];
                               'threshold'     'real'   []              NaN;
                               'unitx'         'string' []              'ms';
                               'chanlocs'      'struct' []              struct('labels', {});
                               'plotsubjects'  'string' { 'on','off' }  'off';
                               'condnames'     'cell'   []              {}; % just for legends
                               'groupnames'    'cell'   []              {}; % just for legends
                               'groupstats'    'cell'   []              {};
                               'condstats'     'cell'   []              {};
                               'interstats'    'cell'   []              {};
                               'titles'        'cell'   []              {};
                               'figure'        'string' { 'on','off' }   'on';
                               'plottopo'      'string' { 'on','off' }   'off';
                               'plotstderr'    'string' { 'on','off','diff','nocurve' }   'off';
                               'plotdiff'      'string' { 'on','off' }   'off';
                               'legend'        { 'string','cell' } { { 'on','off' } {} }  'off';
                               'datatype'      'string' { 'ersp','itc','erp','spec' }    'erp';
                               'plotgroups'    'string' { 'together','apart' }  'apart';
                               'plotmode'      'string' { 'test','condensed' }  'test'; % deprecated
                               'plotconditions'    'string' { 'together','apart' }  'apart' }, 'std_plotcurve');

% opt.figure =  'off'; % test by nima
if isstr(opt), error(opt); end;
opt.singlesubject = 'off';
if length(opt.chanlocs) > 1, opt.plottopo = 'on'; end;
if strcmpi(opt.plottopo, 'on') && size(data{1},3) == 1, opt.singlesubject = 'on'; end;
%if size(data{1},2) == 1,                              opt.singlesubject = 'on'; end;
if all(all(cellfun('size', data, 2)==1))               opt.singlesubject = 'on'; end;
if any(any(cellfun('size', data, 2)==1)), opt.groupstats = {}; opt.condstats = {}; end;
if strcmpi(opt.datatype, 'spec'), opt.unit = 'Hz'; end;
if strcmpi(opt.plotsubjects, 'on')
    opt.plotgroups      = 'apart';
    opt.plotconditions  = 'apart';
end;
if strcmpi(opt.plotconditions, 'together') &&  ~isempty(opt.groupstats), opt.plotconditions = 'apart'; end;
if strcmpi(opt.plotgroups,     'together') &&  ~isempty(opt.condstats) , opt.plotgroups     = 'apart'; end;
if isstr(opt.legend), opt.legend = {}; end;
if isempty(opt.titles), opt.titles = cell(10,10); opt.titles(:) = { '' }; end;
if length(data(:)) == length(opt.legend(:)), 
    opt.legend = reshape(opt.legend, size(data))'; 
    opt.legend(cellfun(@isempty, data)) = []; 
    opt.legend = (opt.legend)';
end;

% color matrix
% -----------------------
onecol  = { 'b' 'b' 'b' 'b' 'b' 'b' 'b' 'b' 'b' 'b' };
manycol = { 'b' 'g' 'm' 'c' 'r' 'k' 'y' 'b' 'g' 'c' 'm' 'r' 'b' 'g' 'c' 'm' 'r' 'b' ...
                   'g' 'c' 'm' 'r' 'b' 'g' 'c' 'm' 'r' 'b' 'g' 'c' 'm' 'r' 'b' 'g' 'c' 'm' };
modifier = { '-' '--' '-.' ':' '-' '--' '-.' ':' '-' '--' '-.' ':'  };
if strcmpi(opt.plotgroups, 'together') || strcmpi(opt.plotconditions, 'together') || strcmpi(opt.figure, 'off')
     col = manycol;
else col = onecol;
end;
nonemptycell = find(~cellfun(@isempty, data));
maxdim = max(length(data(:)), size(data{nonemptycell(1)}, ndims(data{nonemptycell(1)})));
tmpcol = col;
if strcmpi(opt.plotsubjects, 'off')
    
    % both group and conditions together
    if strcmpi(opt.plotconditions, 'together') && strcmpi(opt.plotgroups, 'together')
        dim1 = max(size(data));
        dim2 = min(size(data));
        coldata = col([1:dim1]);
        for iRow = 2:dim2
            coldata(iRow,:) = coldata(1,:);
            for iCol = 1:dim1
                coldata{iRow,iCol} = [ coldata{iRow,iCol} modifier{iRow} ];
            end;
        end;
        if size(coldata,1) ~= size(data,1), coldata = coldata'; end;
    else
        coldata = manycol;
    end;
       
    %coldata = col(mod([0:maxdim-1], length(col))+1);
    %coldata = col(mod([0:maxdim-1], length(col))+1);
    %coldata = reshape(coldata(1:length(data(:))), size(data));
else
    coldata = cell(size(data));
end;

% remove empty entries
% --------------------
datapresent = ~cellfun(@isempty, data);
if size(data,1) > 1, for c = size(data,1):-1:1, if sum(datapresent(c,:)) == 0, data(c,:) = []; coldata(c,:) = []; if ~strcmpi(opt.plotconditions, 'together') opt.titles(c,:) = []; end; if ~isempty(opt.groupstats), opt.groupstats(c) = []; end; end; end; end;
if size(data,2) > 1, for g = size(data,2):-1:1, if sum(datapresent(:,g)) == 0, data(:,g) = []; coldata(:,g) = []; if ~strcmpi(opt.plotgroups    , 'together') opt.titles(:,g) = []; end; if ~isempty(opt.condstats ), opt.condstats( g) = []; end; end; end; end;
if strcmpi(opt.plotsubjects, 'off'), tmpcol = coldata'; tmpcol = tmpcol(:)'; end;

% number of columns and rows to plot
% ----------------------------------
nc = size(data,1);
ng = size(data,2);
if strcmpi(opt.plotgroups, 'together'),      ngplot = 1; else ngplot = ng; end;
if strcmpi(opt.plotconditions,  'together'), ncplot = 1; else ncplot = nc; end;     
if nc >= ng, opt.subplot = 'transpose';
else         opt.subplot = 'normal';
end;
if isempty(opt.condnames)
    for c=1:nc, opt.condnames{c} = sprintf('Cond. %d', c); end;
    if nc == 1, opt.condnames = { '' }; end;
end;
if isempty(opt.groupnames)
    for g=1:ng, opt.groupnames{g} = sprintf('Group. %d', g); end;
    if ng == 1, opt.groupnames = { '' }; end;
end;

% plotting paramters
% ------------------
if ng > 1 && ~isempty(opt.groupstats), addc = 1; else addc = 0; end;
if nc > 1 && ~isempty(opt.condstats ), addr = 1; else addr = 0; end;
if length(opt.threshold) > 1, opt.threshold = opt.threshold(1); end;
if strcmpi(opt.singlesubject, 'off') ...
        && ( ~isempty(opt.condstats) || ~isempty(opt.groupstats) ) % only for curves
    plottag = 0;
    if strcmpi(opt.plotgroups, 'together') && isempty(opt.condstats) && ~isempty(opt.groupstats) && ~isnan(opt.threshold), addc = 0; plottag = 1; end;
    if strcmpi(opt.plotconditions , 'together') && ~isempty(opt.condstats) && isempty(opt.groupstats) && ~isnan(opt.threshold), addr = 0; plottag = 1; end;
    if ~isnan(opt.threshold) && plottag == 0 && strcmpi(opt.figure, 'on')
        disp('Warning: cannot plot condition/group on the same panel while using a fixed');
        disp('         threshold, unless you only compute statistics for ether groups or conditions');
        opt.plotgroups = 'apart';
        opt.plotconditions  = 'apart';
    end;
end;

% resize data to match points x channels x subjects
% or                   points x    1   x components
% -------------------------------------------------
for index = 1:length(data(:))
    if length(opt.chanlocs) ~= size(data{index},2) && (length(opt.chanlocs) == 1 || isempty(opt.chanlocs))
        data{index} = reshape(data{index}, [ size(data{index},1) 1 size(data{index},2) ]);
    end;
end;

% compute significance mask
% --------------------------
if ~isempty(opt.interstats), pinter = opt.interstats{3}; end;

if ~isnan(opt.threshold) && ( ~isempty(opt.groupstats) || ~isempty(opt.condstats) )    
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

% labels
% ------
if strcmpi(opt.unitx, 'ms'), xlab = 'Time (ms)';      ylab = 'Potential (\muV)';
else                         xlab = 'Frequency (Hz)'; ylab = 'Power (10*log_{10}(\muV^{2}))'; 
end;
if ~isnan(opt.threshold), statopt = {  'xlabel' xlab };
else                      statopt = { 'logpval' 'on' 'xlabel' xlab 'ylabel' '-log10(p)' 'ylim' [0 maxplot] };
end;

% adjust figure size
% ------------------
if strcmpi(opt.figure, 'on')
    figure('color', 'w');
    pos = get(gcf, 'position');
    basewinsize = 200/max(nc,ng)*3;
    if strcmpi(opt.plotgroups, 'together') pos(3) = 200*(1+addc);
    else                                   pos(3) = 200*(ng+addc);
    end;
    if strcmpi(opt.plotconditions , 'together') pos(4) = 200*(1+addr);
    else                                        pos(4) = 200*(nc+addr);
    end;
    if strcmpi(opt.subplot, 'transpose'), set(gcf, 'position', [ pos(1) pos(2) pos(4) pos(3)]);
    else                                  set(gcf, 'position', pos);
    end;
else
    opt.subplot = 'noplot';
end;

tmplim = [Inf -Inf];
colcount = 1; % only when plotting all conditions on the same figure
for c = 1:ncplot
    for g = 1:ngplot
        if strcmpi(opt.plotgroups, 'together'),         hdl(c,g)=mysubplot(ncplot+addr, ngplot+addc, 1 + (c-1)*(ngplot+addc), opt.subplot); ci = g;
        elseif strcmpi(opt.plotconditions, 'together'), hdl(c,g)=mysubplot(ncplot+addr, ngplot+addc, g, opt.subplot); ci = c;
        else                                            hdl(c,g)=mysubplot(ncplot+addr, ngplot+addc, g + (c-1)*(ngplot+addc), opt.subplot); ci = 1;
        end;
        
        if ~isempty(data{c,g})

            % read all data from one condition or group
            % -----------------------------------------
            dimreduced_sizediffers = 0;
            if ncplot ~= nc && ngplot ~= ng
                maxdim = max(max(cellfun(@(x)(size(x, ndims(x))), data)));
                for cc = 1:size(data,1)
                    for gg = 1:size(data,2)
                        tmptmpdata = real(data{cc,gg});
                        if cc == 1 && gg == 1, tmpdata = NaN*zeros([size(tmptmpdata,1) size(tmptmpdata,2) maxdim length(data(:))]); end;
                        tmpdata(:,:,1:size(tmptmpdata,3),gg+((cc-1)*ng)) = tmptmpdata;
                    end;
                end;
            elseif ncplot ~= nc % plot conditions together
                for ind = 2:size(data,1), if any(size(data{ind,1}) ~= size(data{1})), dimreduced_sizediffers = 1; end; end;
                for cc = 1:nc
                    tmptmpdata = real(data{cc,g});
                    if dimreduced_sizediffers
                        tmptmpdata = nan_mean(tmptmpdata,ndims(tmptmpdata)); % average across last dim
                    end;
                    if cc == 1 && ndims(tmptmpdata) == 3, tmpdata = zeros([size(tmptmpdata)   nc]); end;
                    if cc == 1 && ndims(tmptmpdata) == 2, tmpdata = zeros([size(tmptmpdata) 1 nc]); end;
                    tmpdata(:,:,:,cc) = tmptmpdata;
                end;
            elseif ngplot ~= ng % plot groups together
                for ind = 2:size(data,2), if any(size(data{1,ind}) ~= size(data{1})), dimreduced_sizediffers = 1; end; end;
                for gg = 1:ng
                    tmptmpdata = real(data{c,gg});
                    if dimreduced_sizediffers
                        tmptmpdata = nan_mean(tmptmpdata,ndims(tmptmpdata));
                    end;
                    if gg == 1 && ndims(tmptmpdata) == 3, tmpdata = zeros([size(tmptmpdata)   ng]); end;
                    if gg == 1 && ndims(tmptmpdata) == 2, tmpdata = zeros([size(tmptmpdata) 1 ng]); end;
                    tmpdata(:,:,:,gg) = tmptmpdata;
                end;
            else
                tmpdata = real(data{c,g});
            end;
            
            % plot difference
            % ---------------
            if ~strcmpi(opt.plotdiff, 'off')
                if ngplot ~= ng || ncplot ~= nc
                    if size(tmpdata,3) == 2
                        tmpdata(:,:,end+1) = tmpdata(:,:,2)-tmpdata(:,:,1);
                        opt.legend{end+1} = [ opt.legend{2} '-' opt.legend{1} ];
                    elseif size(tmpdata,4) == 2
                        tmpdata(:,:,:,end+1) = tmpdata(:,:,:,2)-tmpdata(:,:,:,1);
                        opt.legend{end+1} = [ opt.legend{2} '-' opt.legend{1} ];
                    else
                        disp('Cannot plot difference, more than 2 indep. variable values');
                    end;
                else
                    disp('Cannot plot difference, indep. variable value must be plotted together');
                end;
            end;
            
            if ~isempty(opt.filter), tmpdata = myfilt(tmpdata, 1000/(allx(2)-allx(1)), 0, opt.filter); end;
            
            % plotting options
            % ----------------
            plotopt = { allx };
            % -------------------------------------------------------------
            % tmpdata is of size "points x channels x subject x conditions"
            % or                 "points x   1   x components x conditions"
            % -------------------------------------------------------------
            if ~dimreduced_sizediffers && strcmpi(opt.plotsubjects, 'off') % average accross subjects
                tmpstd = squeeze(real(std(tmpdata,[],3)))/sqrt(size(tmpdata,3)); tmpstd = squeeze(permute(tmpstd, [2 1 3])); tmpdata = squeeze(real(nan_mean(tmpdata,3)));
            end;
            tmpdata = squeeze(permute(tmpdata, [2 1 3 4]));
            % -----------------------------------------------------------------
            % tmpdata is now of size "channels x points x subject x conditions"
            % -----------------------------------------------------------------
            if strcmpi(opt.plottopo, 'on'), highlight = 'background'; else highlight = 'bottom'; end;
            if strcmpi(opt.plotgroups, 'together') &&  isempty(opt.condstats) && ...
                             ~isnan(opt.threshold) && ~isempty(opt.groupstats)
                plotopt = { plotopt{:} 'maskarray' };
                tmpdata = { tmpdata pgroupplot{c}' };
            elseif strcmpi(opt.plotconditions, 'together') &&  isempty(opt.groupstats) && ...
                                     ~isnan(opt.threshold) && ~isempty(opt.condstats)
                plotopt = { plotopt{:} 'maskarray' };
                tmpdata = { tmpdata pcondplot{g}' };
            end;
            plotopt = { plotopt{:} 'highlightmode', highlight };
            if strcmpi(opt.plotsubjects, 'on')
                plotopt = { plotopt{:} 'plotmean' 'on' 'plotindiv' 'on' };
            else
                plotopt = { plotopt{:} 'plotmean' 'off' };
            end;
            plotopt = { plotopt{:} 'ylim' opt.ylim 'xlabel' xlab 'ylabel' ylab };
            if ncplot ~= nc || ngplot ~= ng
                plotopt = { plotopt{:} 'legend' opt.legend };
            end;
            
            if strcmpi(opt.plottopo, 'on') && length(opt.chanlocs) > 1
                metaplottopo(tmpdata, 'chanlocs', opt.chanlocs, 'plotfunc', 'plotcurve', ...
                    'plotargs', { plotopt{:} }, 'datapos', [2 3], 'title', opt.titles{c,g});
            elseif iscell(tmpdata)
                plotcurve( allx, tmpdata{1}, 'colors', tmpcol, 'maskarray', tmpdata{2}, plotopt{3:end}, 'title', opt.titles{c,g});
            else
                if isempty(findstr(opt.plotstderr, 'nocurve'))
                    plotcurve( allx, tmpdata, 'colors', tmpcol, plotopt{2:end}, 'traceinfo', 'on', 'title', opt.titles{c,g});
                end;
                if ~strcmpi(opt.plotstderr, 'off') 
                    if ~dimreduced_sizediffers
                        if ~isempty(findstr(opt.plotstderr, 'diff')), begind = 3; else begind = 1; end;
                        set(gcf, 'renderer', 'OpenGL')
                        for tmpi = begind:size(tmpdata,1)
                            hold on; chandle = fillcurves( allx, tmpdata(tmpi,:)-tmpstd(tmpi,:), tmpdata(tmpi,:)+tmpstd(tmpi,:), tmpcol{tmpi}); hold on;
                            numfaces = size(get(chandle(1), 'Vertices'),1);
                            set(chandle(1), 'FaceVertexCData', repmat([1 1 1], [numfaces 1]), 'Cdatamapping', 'direct', 'facealpha', 0.3, 'edgecolor', 'none');
                        end;
                    else
                        disp('Some conditions have more subjects than others, cannot plot standard error');
                    end;
                end;
            end;
        end;
        
        if strcmpi(opt.plottopo, 'off'), % only non-topographic
            xlim([allx(1) allx(end)]); hold on;
            if isempty(opt.ylim)
                tmp = ylim;
                tmplim = [ min(tmplim(1), tmp(1)) max(tmplim(2), tmp(2)) ];
            else 
                ylim(opt.ylim);
            end;
        end;

        % statistics accross groups
        % -------------------------
        if g == ngplot && ng > 1 && ~isempty(opt.groupstats)            
            if ~strcmpi(opt.plotgroups, 'together') || ~isempty(opt.condstats) || isnan(opt.threshold)
                if strcmpi(opt.plotgroups, 'together'),         mysubplot(ncplot+addr, ngplot+addc, 2 + (c-1)*(ngplot+addc), opt.subplot); ci = g;
                elseif strcmpi(opt.plotconditions, 'together'), mysubplot(ncplot+addr, ngplot+addc, ngplot + 1, opt.subplot); ci = c;
                else                                            mysubplot(ncplot+addr, ngplot+addc, ngplot + 1 + (c-1)*(ngplot+addc), opt.subplot); ci = 1;
                end;
                if strcmpi(opt.plotconditions, 'together'), condnames = 'Conditions'; else condnames = opt.condnames{c}; end;
                if ~isnan(opt.threshold)
                     if strcmpi(opt.plottopo, 'on'), 
                          metaplottopo({zeros(size(pgroupplot{c}')) pgroupplot{c}'}, 'chanlocs', opt.chanlocs, 'plotfunc', 'plotcurve', ...
                              'plotargs', { allx 'maskarray' statopt{:} }, 'datapos', [2 3], 'title', opt.titles{c, g+1});
                     else plotcurve(allx, zeros(size(allx)), 'maskarray', mean(pgroupplot{c},2), 'ylim', [0.1 1], 'title', opt.titles{c, g+1}, statopt{:});
                     end;
                else
                     if strcmpi(opt.plottopo, 'on'), 
                          metaplottopo(pgroupplot{c}', 'chanlocs', opt.chanlocs, 'plotfunc', 'plotcurve', ...
                              'plotargs', { allx statopt{:} }, 'datapos', [2 3], 'title', opt.titles{c, g+1});
                     else plotcurve(allx, mean(pgroupplot{c},2), 'title', opt.titles{c, g+1}, statopt{:});
                     end;
                end;
            end;
        end;
    end;
end;

for g = 1:ng
    % statistics accross conditions
    % -----------------------------
    if ~isempty(opt.condstats) && nc > 1
        if ~strcmpi(opt.plotconditions, 'together') || ~isempty(opt.groupstats) || isnan(opt.threshold)
            if strcmpi(opt.plotgroups, 'together'),         mysubplot(ncplot+addr, ngplot+addc, 1 + c*(ngplot+addc), opt.subplot); ci = g;
            elseif strcmpi(opt.plotconditions, 'together'), mysubplot(ncplot+addr, ngplot+addc, g + ngplot+addc, opt.subplot); ci = c;
            else                                            mysubplot(ncplot+addr, ngplot+addc, g + c*(ngplot+addc), opt.subplot); ci = 1;
            end;
            if strcmpi(opt.plotgroups, 'together'), groupnames = 'Groups'; else groupnames = opt.groupnames{g}; end;
            if ~isnan(opt.threshold)
                 if strcmpi(opt.plottopo, 'on'), 
                      metaplottopo({zeros(size(pcondplot{g}')) pcondplot{g}'}, 'chanlocs', opt.chanlocs, 'plotfunc', 'plotcurve', ...
                          'plotargs', { allx 'maskarray' statopt{:} }, 'datapos', [2 3], 'title', opt.titles{end, g});
                 else plotcurve(allx, zeros(size(allx)), 'maskarray', mean(pcondplot{g},2), 'ylim', [0.1 1], 'title', opt.titles{end, g}, statopt{:});
                 end;
            else
                 if strcmpi(opt.plottopo, 'on'), 
                      metaplottopo(pcondplot{g}', 'chanlocs', opt.chanlocs, 'plotfunc', 'plotcurve', ...
                          'plotargs', { allx statopt{:} }, 'datapos', [2 3], 'title', opt.titles{end, g});
                 else plotcurve(allx, mean(pcondplot{g},2), 'title',  opt.titles{end, g}, statopt{:});
                 end;
            end;
        end;
    end;
end;

% statistics accross group and conditions
% ---------------------------------------
if ~isempty(opt.groupstats) && ~isempty(opt.condstats) && ng > 1 && nc > 1
    mysubplot(ncplot+addr, ngplot+addc, ngplot + 1 + ncplot*(ngplot+addr), opt.subplot);
    if ~isnan(opt.threshold)
         if strcmpi(opt.plottopo, 'on'), 
              metaplottopo({zeros(size(pinterplot')) pinterplot'}, 'chanlocs', opt.chanlocs, 'plotfunc', 'plotcurve', ...
                  'plotargs', { allx 'maskarray' statopt{:} }, 'datapos', [2 3], 'title', opt.titles{end, end});
         else plotcurve(allx, zeros(size(allx)), 'maskarray', mean(pinterplot,2), 'ylim', [0.1 1], 'title', opt.titles{end, end}, statopt{:});
              xlabel(xlab); ylabel('-log10(p)');
        end;
    else
         if strcmpi(opt.plottopo, 'on'), 
              metaplottopo(pinterplot', 'chanlocs', opt.chanlocs, 'plotfunc', 'plotcurve', ...
                  'plotargs', { allx statopt{:} }, 'datapos', [2 3], 'title', opt.titles{end, end});
         else plotcurve(allx, mean(pinterplot,2), 'title', opt.titles{end, end}, statopt{:});
         end;
    end;
end;  

% axis limit
% ----------
for c = 1:ncplot
    for g = 1:ngplot
        if isempty(opt.ylim) && strcmpi(opt.plottopo, 'off')
            set(hdl(c,g), 'ylim', tmplim);
        end;
    end;
end;

if strcmpi(opt.plottopo, 'off') && length(hdl(:)) > 1
    axcopy;
    % remove axis labels (for most but not all)
    % ------------------
    if strcmpi(opt.subplot, 'transpose')
        for c = 1:size(hdl,2)
            for g = 1:size(hdl,1)
                axes(hdl(g,c));
                if c ~= 1 && size(hdl,2) ~=1, xlabel(''); legend off; end;
                if g ~= 1 && size(hdl,1) ~= 1, ylabel(''); legend off; end;
            end;
        end;
    else
        for c = 1:size(hdl,1)
            for g = 1:size(hdl,2)
                axes(hdl(c,g));
                if g ~= 1 && size(hdl,2) ~=1, ylabel(''); legend off; end;
                if c ~= size(hdl,1) && size(hdl,1) ~= 1, xlabel(''); legend off; end;
            end;
        end;
    end;
end;


% mysubplot (allow to transpose if necessary)
% -------------------------------------------
function hdl = mysubplot(nr,nc,ind,subplottype);

    r = ceil(ind/nc);
    c = ind -(r-1)*nc;
    if strcmpi(subplottype, 'transpose'),  hdl = subplot(nc,nr,(c-1)*nr+r);
    elseif strcmpi(subplottype, 'normal'), hdl = subplot(nr,nc,(r-1)*nc+c);
    elseif strcmpi(subplottype, 'noplot'), hdl = gca;
    else error('Unknown subplot type');
    end;

% rapid filtering for ERP
% -----------------------
function tmpdata2 = myfilt(tmpdata, srate, lowpass, highpass)

    tmpdata2 = reshape(tmpdata, size(tmpdata,1), size(tmpdata,2)*size(tmpdata,3)*size(tmpdata,4));
    tmpdata2 = eegfiltfft(tmpdata2',srate, lowpass, highpass)';
    tmpdata2 = reshape(tmpdata2, size(tmpdata,1), size(tmpdata,2), size(tmpdata,3), size(tmpdata,4));
