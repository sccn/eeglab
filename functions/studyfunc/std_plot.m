% std_plot() - plot ERP/spectral traces or ERSP/ITC images a component
%              or channel cluster in a STUDY. Also allows plotting scalp
%              maps.
% Usage:
%          >> std_plot( axvals, data, 'key', 'val', ...)
% Inputs:
%  axvals - [vector or cell array] axis values for the data. 
%           If data are ERSP or ITC images, then this should be a cell 
%           array giving { latencies_in_ms frequencies_in_hz }
%  data  -  [cell array] mean data for each subject group and/or data
%           condition. For example, to plot mean ERPs from a STUDY 
%           for epochs of 800 frames in two conditions from three groups 
%           of 12 subjects:
%
%           >> data = { [800x12] [800x12] [800x12];... % 3 groups, cond 1
%                       [800x12] [800x12] [800x12] };  % 3 groups, cond 2
%           >> std_plot(erp_ms,data);
%
%           By default, parametric statistics are computed across subjects 
%           in the three groups. (group,condition) ERP averages are plotted. 
%           See below and >> help statcond 
%           for more information about the statistical computations.
%
% Optional display parameters:
%  'datatype'    - ['ersp'|'itc'|'erp'|'spec'] data type {default: 'erp'}
%  'channels'    - [cell array] channel names (for titles) {default: all}
%  'condnames'   - [cell array] names of conditions (for titles} 
%                  {default: none}
%  'groupnames'  - [cell array] names of subject groups (for titles)
%                  {default: none}
%  'subject'     - [string] plot subject name (for title)
%
% Statistics options:
%  'groupstats'  - ['on'|'off'] Compute (or not) statistics across groups.
%                  {default: 'off'}
%  'condstats'   - ['on'|'off'] Compute (or not) statistics across groups.
%                  {default: 'off'}
%  'statistics'  - ['param'|'perm'] Type of statistics to use: 'param' for
%                  parametric; 'perm' for permutations {default: 'param'}
%  'naccu'       - [integer] Number of surrogate averges fo accumulate when 
%                  computing permutation-based statistics. For example, to
%                  test p<0.01 use naccu>=200; for p<0.001, use naccu>=2000. 
%                  If a non-NaN 'threshold' is set (see below) and 'naccu' 
%                  is too low, it will be automatically increased. This 
%                  keyword available only from the command line {default:500}
%  'statmode'    - ['subjects'|'trials'] standard statistics are 
%                  'subjects' where the statistics is performed accross
%                  the mean ERSP (or ITC) of single subjects. For 'trials'
%                  statistics, the single-trial data epochs of all subjects
%                  are pooled together. This requires that they were
%                  saved on disk using option 'savetrials', 'on' at the time
%                  of computation. Note that these single-trial data
%                  may use several GB of disk space and that computation 
%                  of 'trials' statistics requires a lot of RAM.
%  'threshold'   - [NaN|real<<1] Significance threshold. NaN -> plot the 
%                  p-values themselves on a different figure. When possible, 
%                  significance regions are indicated below the data.
%                  {default: NaN}
%  'maskdata'    - ['on'|'off'] when threshold is non-NaN and not both 
%                  condition and group statistics are computed, the user 
%                  has the option to mask the data for significance.
%                  {defualt: 'off'}
%
% Curve plotting options (ERP and spectrum):
%  'plotgroups'  - ['together'|'apart'] 'together' -> plot mean results 
%                  for subject groups in the same figure panel in different 
%                  colors. 'apart' -> plot group results on different figure
%                  panels {default: 'apart'}
%  'plotconditions' - ['together'|'apart'] 'together' -> plot mean results 
%                  for data conditions on the same figure panel in different 
%                  colors. 'apart' -> plot conditions on different figure
%                  panel. Note: 'plotgroups' and 'plotconditions' arguments 
%                  cannot both be 'together' {default: 'apart'}
%  'freqrange'   - [min max] spectrum plotting frequency range in Hz. 
%                  {default: all available frequencies}
%  'timerange'   - [min max] ERP plotting latency range in ms. 
%                  {default: all available latencies}
%  'legend'      - ['on'|'off'] turn plot legend on/off {default: 'off'}
%  'plotmode'    - ['normal'|'condensed'] statistics plotting mode:
%                  'condensed' -> plot statistics under the curves 
%                  (when possible); 'normal' -> plot them in separate 
%                  axes {default: 'normal'}
% 'plotsubjects' - ['on'|'off'] overplot traces for individual components
%                  or channels {default: 'off'}
%  'plottopo'    - [real] plot a scalp map for the specified latency (ms)
%                  or frequency (Hz) {default: no topoplot}
%  'ylim'        - [min max] ordinate limits for ERP and spectrum plots
%                  {default: all available data}
%
% ITC/ERSP image plotting options:
%  'plottopo'    - [latency frequency] plot a scalp map at specific 
%                  latency (in ms) and frequency (in Hz).
%  'tftopoopt'   - [cell array] tftopo() plotting options (ERSP and ITC)
%  'caxis'       - [min max] color axis (ERSP, ITC, scalp maps)
%
% Scalp map plotting options:
%  'chanlocs'    - [struct] channel location structure
%
% Author: Arnaud Delorme, CERCO, CNRS, 2006-
%
% See also: pop_erspparams(), pop_erpparams(), pop_specparams(), statcond()

% $Log: not supported by cvs2svn $
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

function [pgroup, pcond, pinter] = std_plot(allx, data, varargin)

pgroup = [];
pcond  = [];
pinter = [];
if nargin < 2
    help std_plot;
    return;
end;

opt = finputcheck( varargin, { 'channels'    'cell'   []              {};
                               'caxis'       'real'   []              [];
                               'ylim'        'real'   []              [];
                               'condnames'    'cell'   []              {};
                               'groupnames'   'cell'   []              {};
                               'tftopoopt'   'cell'   []              {};
                               'threshold'   'real'   []              NaN;
                               'plottopo'       'real'   []              [];
                               'naccu'       'integer' []             500;
                               'unitx'       'string' []              'ms'; % just for titles
                               'subject'     'string' []              '';   % just for titles
                               'chanlocs'    'struct' []              struct('labels', {});
                               'plotsubjects' 'string' { 'on' 'off' }  'off';
                               'groupstats'   'string' { 'on' 'off' }   'off';
                               'condstats'    'string' { 'on' 'off' }   'off';
                               'maskdata'    'string' { 'on' 'off' }   'off';
                               'legend'      'string' { 'on' 'off' }   'off';
                               'timerange'   'real'   []               [];
                               'datatype'    'string' { 'ersp' 'itc' 'erp' 'spec' }    'erp';
                               'plotgroups'   'string' { 'together' 'apart' }  'apart';
                               'plotconditions'    'string' { 'together' 'apart' }  'apart';
                               'plotmode'    'string' { 'normal' 'condensed' }  'normal';
                               'statistics'  'string' { 'param' 'perm' }       'param';
                               'statmode'    'string' { 'subjects' 'common' 'trials' } 'subjects'}, 'std_erpmaskdata');
                           
if isstr(opt), error(opt); end;
if strcmpi(opt.datatype, 'spec'), g.unit = 'Hz'; end;
if strcmpi(opt.plotsubjects, 'on')
    opt.plotgroups = 'apart';
    opt.plotconditions  = 'apart';
end;
onecol  = { 'b' 'b' 'b' 'b' 'b' 'b' 'b' 'b' 'b' 'b' };
manycol = { 'b' 'r' 'g' 'k' 'c' 'y' };
                            
nc = size(data,1);
ng = size(data,2);
if nc >= ng, opt.transpose = 'on';
else         opt.transpose = 'off';
end;
if isempty(opt.condnames)
    for c=1:nc, opt.condnames{c} = sprintf('Cond. %d', c); end;
    if nc == 1, opt.condnames = { '' }; end;
end;
if isempty(opt.groupnames)
    for g=1:ng, opt.groupnames{g} = sprintf('Group. %d', g); end;
    if ng == 1, opt.groupnames = { '' }; end;
end;

% condensed plot
% --------------
if strcmpi(opt.plotmode, 'condensed') 
    if strcmpi(opt.datatype, 'erp') | strcmpi(opt.datatype, 'spec')
        leg = {};
        for c = 1:nc
            for g = 1:ng
                plot( allx, real(mean(mean(data{c,g},2),3)), manycol{(c-1)*ng+g});
                xlim([allx(1) allx(end)]); hold on;
                if ~isempty(opt.ylim), ylim(opt.ylim); end;
                if isempty(opt.condnames{c}) | isempty(opt.groupnames{g}), leg{(c-1)*ng+g} = [ opt.condnames{c} opt.groupnames{g} ];
                else                                                     leg{(c-1)*ng+g} = [ opt.condnames{c} ', ' opt.groupnames{g} ];
                end;
            end;
        end;
        if strcmpi(opt.legend, 'on'), legend(leg{:}); end;
    else
        meanplot = zeros(size(data{1},1), size(data{1},2));
        for c = 1:nc
            for g = 1:ng
                meanplot = meanplot + mean(data{c,g},3)/nc/ng;
            end;
        end;
        options = { 'chanlocs', opt.chanlocs, 'electrodes', 'off', 'cbar', 'off', ...
                'cmode', 'separate', 'logfreq', 'native', opt.tftopoopt{:} };
        tftopo( meanplot', allx{2}, allx{1}, 'title', [ 'Mean ' upper(opt.datatype) ' for all group/cond' ], options{:}); 
    end;
    return;
end;

% plotting paramters
% ------------------
if ng > 1 & strcmpi(opt.groupstats, 'on'), addc = 1; else addc = 0; end;
if nc > 1 & strcmpi(opt.condstats , 'on'), addr = 1; else addr = 0; end;
if isempty(opt.plottopo) % only for curves
    plottag = 0;
    if strcmpi(opt.plotgroups, 'together') & strcmpi(opt.condstats, 'off') & strcmpi(opt.groupstats, 'on' ) & ~isnan(opt.threshold), addc = 0; plottag = 1; end;
    if strcmpi(opt.plotconditions , 'together') & strcmpi(opt.condstats, 'on' ) & strcmpi(opt.groupstats, 'off') & ~isnan(opt.threshold), addr = 0; plottag = 1; end;
    if ~isnan(opt.threshold) & plottag == 0
        disp('Warning: cannot plot condition/group on the same panel while using a fixed');
        disp('         threshold, unless you only compute statistics for ether groups or conditions');
        opt.plotgroups = 'apart';
        opt.plotconditions  = 'apart';
    end;
end;

if ~isempty(opt.plottopo)
    % ----------------------------
    % plot scalp maps for baseline    
    % ----------------------------
    if length(opt.plottopo) < 1, opt.plottopo(2) = opt.plottopo(1); end;
    [tmp ti1] = min(abs(allx-opt.plottopo(1)));
    [tmp ti2] = min(abs(allx-opt.plottopo(2)));
    for index = 1:length(data(:))
        data{index} = mean(data{index}(ti1:ti2,:,:),1);
    end;
end;

% compute significance mask
% --------------------------
if strcmpi(opt.condstats, 'on') & nc > 1
    for g = 1:ng
        [F df pval] = condstat(data(:,g), 'mode', opt.statistics, 'naccu', opt.naccu); 
        pcond{g} = squeeze(pval);
    end;
else
    pcond = {};
end;
if strcmpi(opt.groupstats, 'on') & ng > 1
    for c = 1:nc
        [F df pval] = condstat(data(c,:), 'mode', opt.statistics, 'naccu', opt.naccu); 
        pgroup{c} = squeeze(pval);
    end;
else
    pgroup = {};
end;
if ( strcmpi(opt.groupstats, 'on') | strcmpi(opt.condstats, 'on') ) & ng > 1 & nc > 1
    [F df pval] = condstat(data, 'mode', opt.statistics, 'naccu', opt.naccu);
    pinter      = squeeze(pval{3});
else
    pinter = [];
end;
if ~isnan(opt.threshold)    
    % applying threshold
    % ------------------
    for ind = 1:length(pcond),  pcondplot{ind}  = pcond{ind}  < opt.threshold; end;
    for ind = 1:length(pgroup), pgroupplot{ind} = pgroup{ind} < opt.threshold; end;
    pinterplot = pinter < opt.threshold;
    maxplot = 1;
else
    warning off;
    for ind = 1:length(pcond),  pcondplot{ind}  = -log10(pcond{ind}); end;
    for ind = 1:length(pgroup), pgroupplot{ind} = -log10(pgroup{ind}); end;
    pinterplot = -log10(pinter);
    maxplot = 3;
    warning on;
end;

% plotting all conditions
% -----------------------
if isempty(opt.plottopo) & ( strcmpi(opt.datatype, 'erp') | strcmpi(opt.datatype, 'spec') )

    if strcmpi(opt.plotgroups, 'together') & strcmpi(opt.plotconditions,  'together')
        error('Cannot plot both conditions and groups on the same plot');
    end;
    if strcmpi(opt.plotgroups, 'together'), ngplot = 1; else ngplot = ng; end; 
    if strcmpi(opt.plotconditions,  'together'), ncplot = 1; else ncplot = nc; end;     
    if strcmpi(opt.plotgroups, 'together') | strcmpi(opt.plotconditions, 'together')
         col = manycol;
         leg = 'on';
    else
         col = onecol;
         leg = 'off';
    end;

    % adjust figure size
    % ------------------
    figure('color', 'w');
    pos = get(gcf, 'position');
    basewinsize = 200/max(nc,ng)*3;
    if strcmpi(opt.plotgroups, 'together') pos(3) = 200*(1+addc);
    else                                  pos(3) = 200*(ng+addc);
    end;
    if strcmpi(opt.plotconditions , 'together') pos(4) = 200*(1+addr);
    else                                  pos(4) = 200*(nc+addr);
    end;
    if strcmpi(opt.transpose, 'on'), set(gcf, 'position', [ pos(1) pos(2) pos(4) pos(3)]);
    else                             set(gcf, 'position', pos);
    end;

    tmplim = [Inf -Inf];
    for c = 1:nc
        for g = 1:ng
            if strcmpi(opt.plotgroups, 'together'),    hdl(c,g)=mysubplot(ncplot+addr, ngplot+addc, 1 + (c-1)*(ngplot+addc), opt.transpose); ci = g;
            elseif strcmpi(opt.plotconditions, 'together'), hdl(c,g)=mysubplot(ncplot+addr, ngplot+addc, g, opt.transpose); ci = c;
            else                                      hdl(c,g)=mysubplot(ncplot+addr, ngplot+addc, g + (c-1)*(ngplot+addc), opt.transpose); ci = 1;
            end;
            if strcmpi(opt.plotsubjects, 'on')
                tmp = plot( allx, real(data{c,g}), 'k'); set(tmp, 'color', [0.5 0.5 0.5]); hold on;
                plot( allx, real(mean(mean(data{c,g},2),3)), 'k', 'linewidth', 2);
            else
                plot( allx, real(mean(mean(data{c,g},2),3)), col{ci});
            end;
            xlim([allx(1) allx(end)]); hold on;
            if isempty(opt.ylim)
                tmp = ylim;
                tmplim = [ min(tmplim(1), tmp(1)) max(tmplim(2), tmp(2)) ];
            else 
                ylim(opt.ylim);
            end;

            % statistics accross groups
            % -------------------------
            if g == ng & ng > 1 & strcmpi(opt.groupstats, 'on')
                if strcmpi(opt.plotgroups, 'together') & strcmpi(opt.condstats, 'off') & ~isnan(opt.threshold)
                     plotcurve(allx, zeros(size(allx)), 'maskarray', mean(pgroupplot{c},2), 'highlightmode', 'bottom');
                else
                    if strcmpi(opt.plotgroups, 'together'),    mysubplot(ncplot+addr, ngplot+addc, 2 + (c-1)*(ngplot+addc), opt.transpose); ci = g;
                    elseif strcmpi(opt.plotconditions, 'together'), mysubplot(ncplot+addr, ngplot+addc, ngplot + 1, opt.transpose); ci = c;
                    else                                      mysubplot(ncplot+addr, ngplot+addc, ngplot + 1 + (c-1)*(ngplot+addc), opt.transpose); ci = 1;
                    end;
                    if strcmpi(opt.plotconditions, 'together'), condnames = 'Conditions'; else condnames = opt.condnames{c}; end;
                    if ~isnan(opt.threshold)
                         plotcurve(allx, zeros(size(allx)), 'maskarray', mean(pgroupplot{c},2));
                         xlim([allx(1) allx(end)]); ylim([0 maxplot]); set(gca, 'ytick', []); title(sprintf('%s (p<%.4f)', condnames, opt.threshold));
                    else plot( allx, mean(pgroupplot{c},2), col{ci});  hold on;
                         xlim([allx(1) allx(end)]); ylim([0 maxplot]); title([ condnames ' (p-value)' ]);
                         set(gca, 'yticklabel', round(10.^-get(gca, 'ytick')*1000)/1000, 'ydir', 'reverse', 'ytickmode', 'manual');
                    end;
                end;
            end;
        end;

    end;

    for g = 1:ng
        % statistics accross conditions
        % -----------------------------
        if strcmpi(opt.condstats, 'on') & nc > 1
            if strcmpi(opt.plotconditions, 'together') & strcmpi(opt.groupstats, 'off') & ~isnan(opt.threshold)
                 axes(hdl(c,g));
                 plotcurve(allx, zeros(size(allx)), 'maskarray', mean(pcondplot{g},2), 'highlightmode', 'bottom');
                 title(opt.groupnames{g}); 
            else
                if strcmpi(opt.plotgroups, 'together'),    mysubplot(ncplot+addr, ngplot+addc, 1 + c*(ngplot+addc), opt.transpose); ci = g;
                elseif strcmpi(opt.plotconditions, 'together'), mysubplot(ncplot+addr, ngplot+addc, g + ngplot+addc, opt.transpose); ci = c;
                else                                      mysubplot(ncplot+addr, ngplot+addc, g + c*(ngplot+addc), opt.transpose); ci = 1;
                end;
                if strcmpi(opt.plotgroups, 'together'), groupnames = 'Groups'; else groupnames = opt.groupnames{g}; end;
                if ~isnan(opt.threshold)
                     plotcurve(allx, zeros(size(allx)), 'maskarray', mean(pcondplot{g},2));
                     xlim([allx(1) allx(end)]); set(gca, 'ytick', []); ylim([0 maxplot]); title(sprintf('%s (p<%.4f)', groupnames, opt.threshold));
                else plot( allx, mean(pcondplot{g},2), col{ci}); hold on;
                     xlim([allx(1) allx(end)]); ylim([0 maxplot]); title([ groupnames ' (p-value)' ]);
                     set(gca, 'yticklabel', round(10.^-get(gca, 'ytick')*1000)/1000, 'ydir', 'reverse');
                end;
            end;
        end;
    end;

    % axis limit and legend
    % ---------------------
    for c = 1:nc
        for g = 1:ng
            axes(hdl(c,g));
            if isempty(opt.ylim)
                ylim(tmplim);
            end;
            if strcmpi(opt.plotgroups, 'together'),                       fig_title = opt.condnames{c};
            elseif strcmpi(opt.plotconditions, 'together'),                    fig_title = opt.groupnames{g};
            elseif isempty(opt.condnames{c}) | isempty(opt.groupnames{g}), fig_title = [ opt.condnames{c} opt.groupnames{g} ];
            else                                                         fig_title = [ opt.condnames{c} ', ' opt.groupnames{g} ];
            end;
            if ~isempty(opt.subject), fig_title = [ opt.subject ', ' fig_title ];
            end;
            title(fig_title); 
            if strcmpi(leg, 'on') & c == nc & g == ng 
                if strcmpi(opt.plotgroups, 'together'), legend( opt.groupnames );
                else                                   legend( opt.condnames );
                end;
            end;
        end;
    end;

    % statistics accross group and conditions
    % ---------------------------------------
    if strcmpi(opt.groupstats, 'on') & strcmpi(opt.condstats, 'on') & ng > 1 & nc > 1
        mysubplot(ncplot+addr, ngplot+addc, ngplot + 1 + ncplot*(ngplot+addr), opt.transpose);
        if ~isnan(opt.threshold)
             plotcurve(allx, zeros(size(allx)), 'maskarray', mean(pinterplot,2));
             xlim([allx(1) allx(end)]); set(gca, 'ytick', []); ylim([0 maxplot]); title(sprintf('Interaction (p<%.4f)', opt.threshold));
        else plot( allx, mean(pinterplot,2), col{ci}); 
             xlim([allx(1) allx(end)]); ylim([0 maxplot]); title('Interaction (p-value)');
             set(gca, 'yticklabel', round(10.^-get(gca, 'ytick')*1000)/1000, 'ydir', 'reverse');
        end;
    end;  
    axcopy;
elseif isempty(opt.plottopo)
    
    % -------------------------------
    % masking for significance of not
    % -------------------------------
    statmask = 0;
    if strcmpi(opt.maskdata, 'on') & ~isnan(opt.threshold) & ...
            (strcmpi(opt.condstats, 'off') | strcmpi(opt.condstats, 'on'))
        addc = 0; addr = 0; statmask = 1;
    end;
        
    % -------------------------
    % plot time/frequency image
    % -------------------------
    options = { 'chanlocs', opt.chanlocs, 'electrodes', 'off', 'cbar', 'off', ...
                'cmode', 'separate', 'logfreq', 'native', opt.tftopoopt{:} };
    figure('color', 'w');
    tmpc = single([inf -inf]);
    for c = 1:nc
        for g = 1:ng
            hdl(c,g) = mysubplot(nc+addr, ng+addc, g + (c-1)*(ng+addc), opt.transpose);
            fig_title = [ opt.condnames{c} ', ' opt.groupnames{g} ];
            tmpplot = mean(data{c,g},3);
            if statmask, 
                if strcmpi(opt.condstats, 'on'), tmpplot(find(pcondplot{g}(:) == 0)) = 0;
                else                            tmpplot(find(pgroupplot{c}(:) == 0)) = 0;
                end;
            end;
            tftopo( tmpplot', allx{2}, allx{1}, 'title', fig_title, options{:}); 
            if isempty(opt.caxis)
                tmpc = [ min(min(tmpplot(:)), tmpc(1)) max(max(tmpplot(:)), tmpc(2)) ];
            else 
                caxis(opt.caxis);
            end;

            % statistics accross groups
            % -------------------------
            if g == ng && ng > 1 & strcmpi(opt.groupstats, 'on') & ~statmask
                hdl(c,g+1) = mysubplot(nc+addr, ng+addc, g + 1 + (c-1)*(ng+addc), opt.transpose);
                if isnan(opt.threshold), tmp_title = sprintf('%s (p-value)', opt.condnames{c});
                else                     tmp_title = sprintf('%s (p<%.4f)',  opt.condnames{c}, opt.threshold);
                end;
                tftopo( pgroupplot{c}', allx{2}, allx{1}, 'title', tmp_title, options{:});
                caxis([-maxplot maxplot]);
            end;
            
        end;

    end;
    for g = 1:ng
        % statistics accross conditions
        % -----------------------------
        if strcmpi(opt.condstats, 'on') & ~statmask && nc > 1
            hdl(nc+1,g) = mysubplot(nc+addr, ng+addc, g + c*(ng+addc), opt.transpose);
            if isnan(opt.threshold), tmp_title = sprintf('%s (p-value)', opt.groupnames{g});
            else                     tmp_title = sprintf('%s (p<%.4f)',  opt.groupnames{g}, opt.threshold);
            end;
            tftopo( pcondplot{g}', allx{2}, allx{1}, 'title', tmp_title, options{:});
            caxis([-maxplot maxplot]);
        end;
    end;
    
    % color scale
    % -----------
    if isempty(opt.caxis)
        tmpc = [-max(abs(tmpc)) max(abs(tmpc))];
        for c = 1:nc
            for g = 1:ng
                axes(hdl(c,g));
                caxis(tmpc);
            end;
        end;
    end;
    
    % statistics accross group and conditions
    % ---------------------------------------
    if strcmpi(opt.groupstats, 'on') & strcmpi(opt.condstats, 'on') && ng > 1 && nc > 1
        hdl(nc+1,ng+1) = mysubplot(nc+addr, ng+addc, g + 1 + c*(ng+addr), opt.transpose);
        if isnan(opt.threshold), tmp_title = 'Interaction (p-value)';
        else                     tmp_title = sprintf('Interaction (p<%.4f)', opt.threshold);
        end;
        tftopo( pinterplot', allx{2}, allx{1}, 'title', tmp_title, options{:});
        caxis([-maxplot maxplot]);
        ylabel('');
    end;    
    
    % color bars
    % ----------
    axes(hdl(nc,ng)); 
    cbar_standard(opt.datatype, ng);
    if nc ~= size(hdl,1) | ng ~= size(hdl,2)
        axes(hdl(end,end));
        cbar_signif(ng, maxplot);
    end;
else    
    
    % topoplot
    % --------
    figure('color', 'w');
    tmpc = [inf -inf];
    for c = 1:nc
        for g = 1:ng
            hdl(c,g) = mysubplot(nc+addr, ng+addc, g + (c-1)*(ng+addc), opt.transpose);
            if opt.plottopo(1) == opt.plottopo(2), fig_title = [ opt.condnames{c} ', ' opt.groupnames{g} ', ' num2str(opt.plottopo(1)) opt.unitx];
            else                             fig_title = [ opt.condnames{c} ', ' opt.groupnames{g} ', ' num2str(opt.plottopo(1)) '-' num2str(opt.plottopo(2)) opt.unitx];
            end;
            tmpplot = double(mean(data{c,g},3));
            topoplot( tmpplot, opt.chanlocs);
            title(fig_title); 
            if isempty(opt.caxis)
                tmpc = [ min(min(tmpplot), tmpc(1)) max(max(tmpplot), tmpc(2)) ];
            else 
                caxis(opt.caxis);
            end;

            % statistics accross groups
            % -------------------------
            if g == ng & ng > 1 & strcmpi(opt.groupstats, 'on')
                hdl(c,g+1) = mysubplot(nc+addr, ng+addc, g + 1 + (c-1)*(ng+addc), opt.transpose);
                topoplot( pgroupplot{c}, opt.chanlocs);
                if isnan(opt.threshold), title(sprintf('%s (p-value)', opt.condnames{c}));
                else                     title(sprintf('%s (p<%.4f)',  opt.condnames{c}, opt.threshold));
                end;
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
        if strcmpi(opt.condstats, 'on') & nc > 1
            hdl(nc+1,g) = mysubplot(nc+addr, ng+addc, g + c*(ng+addc), opt.transpose);
            topoplot( pcondplot{g}, opt.chanlocs);
            if isnan(opt.threshold), title(sprintf('%s (p-value)', opt.groupnames{g}));
            else                     title(sprintf('%s (p<%.4f)',  opt.groupnames{g}, opt.threshold));
            end;
            caxis([-maxplot maxplot]);
        end;
    end;

    % statistics accross group and conditions
    % ---------------------------------------
    if strcmpi(opt.groupstats, 'on') & strcmpi(opt.condstats, 'on') & ng > 1 & nc > 1
        hdl(nc+1,ng+1) = mysubplot(nc+addr, ng+addc, g + 1 + c*(ng+addr), opt.transpose);
        topoplot( pinterplot, opt.chanlocs);
        if isnan(opt.threshold), title('Interaction (p-value)');
        else                     title(sprintf('Interaction (p<%.4f)', opt.threshold));
        end;
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

end;

% remove axis labels
% ------------------
for c = 1:size(hdl,1)
    for g = 1:size(hdl,2)
        if g ~= 1 & size(hdl,2) ~=1, ylabel(''); end;
        if c ~= size(hdl,1) & size(hdl,1) ~= 1, xlabel(''); end;
    end;
end;

% mysubplot (allow to transpose if necessary)
% -------------------------------------------
function hdl = mysubplot(nr,nc,ind,transp);

    r = ceil(ind/nc);
    c = ind -(r-1)*nc;
    if strcmpi(transp, 'on'), hdl = subplot(nc,nr,(c-1)*nr+r);
    else                      hdl = subplot(nr,nc,(r-1)*nc+c);
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
