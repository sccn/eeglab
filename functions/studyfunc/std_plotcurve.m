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
%           for more information about the statistical computations.
%
% Optional display parameters:
%  'datatype'    - ['erp'|'spec'] data type {default: 'erp'}
%  'channels'    - [cell array] channel names (for titles) {default: all}
%  'condnames'   - [cell array] names of conditions (for titles} 
%                  {default: none}
%  'groupnames'  - [cell array] names of subject groups (for titles)
%                  {default: none}
%  'subject'     - [string] plot subject name (for title)
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
%  'plotmode'    - ['normal'|'condensed'] statistics plotting mode:
%                  'condensed' -> plot statistics under the curves 
%                  (when possible); 'normal' -> plot them in separate 
%                  axes {default: 'normal'}
% 'plotsubjects' - ['on'|'off'] overplot traces for individual components
%                  or channels {default: 'off'}
% 'ylim'         - [min max] ordinate limits for ERP and spectrum plots
%                  {default: all available data}
%
% Scalp map plotting options:
%  'chanlocs'    - [struct] channel locations structure
%
% Author: Arnaud Delorme, CERCO, CNRS, 2006-
%
% See also: pop_erspparams(), pop_erpparams(), pop_specparams(), statcond()

% $Log: not supported by cvs2svn $
% Revision 1.15  2008/04/16 18:03:41  arno
% use nan_mean in case some subjects are missing data
%
% Revision 1.14  2007/09/11 10:49:57  arno
% fix numerous small display problems
% and crash
%
% Revision 1.13  2007/08/25 02:03:26  scott
% spelling
%
% Revision 1.12  2007/08/25 01:05:03  arno
% nothing
%
% Revision 1.11  2007/08/14 19:29:55  nima
% _
%
% Revision 1.10  2007/08/09 18:29:03  arno
% plotting statistics
%
% Revision 1.9  2007/08/03 23:34:23  arno
% fix labels
%
% Revision 1.8  2007/04/06 01:58:05  arno
% legend business
%
% Revision 1.7  2007/04/06 01:47:23  arno
% better axis label
%
% Revision 1.6  2007/04/06 01:28:35  arno
% same
%
% Revision 1.5  2007/04/06 01:27:05  arno
% fixed the axis labels
%
% Revision 1.4  2007/03/14 00:56:17  arno
% plotting ERP for multiple components
%
% Revision 1.3  2007/02/28 12:04:38  arno
% do not output statistics
%
% Revision 1.1  2007/01/26 18:08:28  arno
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

function std_plotcurve(allx, data, varargin)

pgroup = [];
pcond  = [];
pinter = [];
if nargin < 2
    help std_plotcurve;
    return;
end;

opt = finputcheck( varargin, { 'channels'    'cell'   []              {};
                               'ylim'        'real'   []              [];
                               'filter'      'real'   []              [];
                               'condnames'   'cell'   []              {};
                               'groupnames'  'cell'   []              {};
                               'compinds'    'cell'   []              {};
                               'threshold'   'real'   []              NaN;
                               'topovals'       'real'   []              []; % same as above
                               'naccu'       'integer' []             500;
                               'unitx'       'string' []              'ms'; % just for titles
                               'subject'     'string' []              '';   % just for titles
                               'chanlocs'    'struct' []              struct('labels', {});
                               'plotsubjects' 'string' { 'on' 'off' }  'off';
                               'groupstats'   'cell'   []              {};
                               'condstats'    'cell'   []              {};
                               'interstats'   'cell'   []              {};
                               'plottopo'     'string' { 'on' 'off' }   'off';
                               'figure'       'string' { 'on' 'off' }   'on';
                               'legend'      'string' { 'on' 'off' }   'off';
                               'datatype'    'string' { 'ersp' 'itc' 'erp' 'spec' }    'erp';
                               'plotgroups'   'string' { 'together' 'apart' }  'apart';
                               'plotconditions'    'string' { 'together' 'apart' }  'apart';
                               'plotmode'    'string' { 'normal' 'condensed' }  'normal';
                               'statistics'  'string' { 'param' 'perm' }       'param';
                               'statmode'    'string' { 'subjects' 'common' 'trials' } 'subjects'}, 'std_erpmaskdata');

% opt.figure =  'off'; % test by nima
% opt.plotmode =  'condensed';
if isstr(opt), error(opt); end;
opt.singlesubject = 'off';
if strcmpi(opt.plottopo, 'on') & size(data{1},3) == 1, opt.singlesubject = 'on'; end;
%if size(data{1},2) == 1,                              opt.singlesubject = 'on'; end;
if all(all(cellfun('size', data, 2)==1))               opt.singlesubject = 'on'; end;
if any(any(cellfun('size', data, 2)==1)), opt.groupstats = {}; opt.condstats = {}; end;
if ~isempty(opt.compinds), if length(opt.compinds{1}) > 1, opt.compinds = {}; end; end;
if strcmpi(opt.datatype, 'spec'), opt.unit = 'Hz'; end;
if strcmpi(opt.plotsubjects, 'on')
    opt.plotgroups      = 'apart';
    opt.plotconditions  = 'apart';
end;
onecol  = { 'b' 'b' 'b' 'b' 'b' 'b' 'b' 'b' 'b' 'b' };
manycol = { 'b' 'r' 'g' 'k' 'c' 'y' };

nc = size(data,1);
ng = size(data,2);
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

% condensed plot (deprecated)
% --------------
if strcmpi(opt.plotmode, 'condensed') 
    opt.plotgroups     = 'together';
    opt.plotconditions = 'together';
    if ~isempty(opt.condstats) & ~isempty(opt.groupstats)
        opt.condstats      = {};
        opt.groupstats     = {};
    end;
end;

% plotting paramters
% ------------------
if ng > 1 & ~isempty(opt.groupstats), addc = 1; else addc = 0; end;
if nc > 1 & ~isempty(opt.condstats ), addr = 1; else addr = 0; end;
if isempty(opt.topovals) & strcmpi(opt.singlesubject, 'off') & ~strcmpi(opt.plotmode, 'condensed') ...
        & ( ~isempty(opt.condstats) | ~isempty(opt.groupstats) ) % only for curves
    plottag = 0;
    if strcmpi(opt.plotgroups, 'together') & isempty(opt.condstats) & ~isempty(opt.groupstats) & ~isnan(opt.threshold), addc = 0; plottag = 1; end;
    if strcmpi(opt.plotconditions , 'together') & ~isempty(opt.condstats) & isempty(opt.groupstats) & ~isnan(opt.threshold), addr = 0; plottag = 1; end;
    if ~isnan(opt.threshold) & plottag == 0 & strcmpi(opt.figure, 'on')
        disp('Warning: cannot plot condition/group on the same panel while using a fixed');
        disp('         threshold, unless you only compute statistics for ether groups or conditions');
        opt.plotgroups = 'apart';
        opt.plotconditions  = 'apart';
    end;
end;

% compute significance mask
% --------------------------
if ~isempty(opt.interstats), pinter = opt.interstats{3}; end;

if ~isnan(opt.threshold)    
    % applying threshold
    % ------------------
    for ind = 1:length(opt.condstats),  pcondplot{ind}  = opt.condstats{ind}  < opt.threshold; end;
    for ind = 1:length(opt.groupstats), pgroupplot{ind} = opt.groupstats{ind} < opt.threshold; end;
    if ~isempty(pinter), pinterplot = pinter < opt.threshold; end;
    maxplot = 1;
else
    warning off;
    for ind = 1:length(opt.condstats),  pcondplot{ind}  = -log10(opt.condstats{ind}); end;
    for ind = 1:length(opt.groupstats), pgroupplot{ind} = -log10(opt.groupstats{ind}); end;
    if ~isempty(pinter), pinterplot = -log10(pinter); end;
    maxplot = 3;
    warning on;
end;

% plotting all conditions
% -----------------------
if strcmpi(opt.plotgroups, 'together'),      ngplot = 1; else ngplot = ng; end; 
if strcmpi(opt.plotconditions,  'together'), ncplot = 1; else ncplot = nc; end;     
if strcmpi(opt.plotgroups, 'together') | strcmpi(opt.plotconditions, 'together') | strcmpi(opt.figure, 'off')
     col = manycol;
     leg = 'on';
else
     col = onecol;
     leg = 'off';
end;

% labels
% ------
if strcmpi(opt.unitx, 'ms'), xlab = 'Time (ms)';      ylab = 'Potential (\muV)';
else                         xlab = 'Frequency (Hz)'; ylab = 'Power (10*log_{10}(\muV^{2}/Hz))'; 
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
for c = 1:ncplot
    for g = 1:ngplot
        if strcmpi(opt.plotgroups, 'together'),         hdl(c,g)=mysubplot(ncplot+addr, ngplot+addc, 1 + (c-1)*(ngplot+addc), opt.subplot); ci = g;
        elseif strcmpi(opt.plotconditions, 'together'), hdl(c,g)=mysubplot(ncplot+addr, ngplot+addc, g, opt.subplot); ci = c;
        else                                            hdl(c,g)=mysubplot(ncplot+addr, ngplot+addc, g + (c-1)*(ngplot+addc), opt.subplot); ci = 1;
        end;
        
        if ~isempty(data{c,g})

            % read all data from one condition or group
            % -----------------------------------------
            if ncplot ~= nc & ngplot ~= ng
                for cc = 1:length(data,1)
                    for gg = 1:length(data,2)
                        tmptmpdata = real(data{cc,gg});
                        if cc == 1, tmpdata = zeros([size(tmptmpdata) length(data(:))]); end;
                        if ndims(tmptmpdata) == 3, tmpdata(:,:,:,gg+((c-1)*ng)) = tmptmpdata;
                        else                       tmpdata(:,:,gg+((cc-1)*ng))  = tmptmpdata;
                        end;
                    end;
                    if isempty(opt.condnames{c}) | isempty(opt.groupnames{g}), leg{(cc-1)*ng+gg} = [ opt.condnames{c} opt.groupnames{g} ];
                    else                                                       leg{(cc-1)*ng+gg} = [ opt.condnames{c} ', ' opt.groupnames{g} ];
                    end;
                end;
            elseif ncplot ~= nc
                for cc = 1:nc
                    tmptmpdata = real(data{cc,g});
                    if cc == 1, tmpdata = zeros([size(tmptmpdata) nc]); end;
                    if ndims(tmptmpdata) == 3, tmpdata(:,:,:,cc) = tmptmpdata;
                    else                       tmpdata(:,:,cc)   = tmptmpdata;
                    end;
                end;
                leg = opt.condnames;
            elseif ngplot ~= ng
                samesize = 1;
                for ind = 2:size(data,2), if any(size(data{1,ind}) ~= size(data{1})), samesize = 0; end; end;
                for gg = 1:ng
                    tmptmpdata = real(data{c,gg});
                    if ~samesize
                        tmptmpdata = nan_mean(tmptmpdata,2);
                    end;
                    if gg == 1, tmpdata = zeros([size(tmptmpdata) ng]); end;
                    if ndims(tmptmpdata) == 3, tmpdata(:,:,:,gg) = tmptmpdata;
                    else                       tmpdata(:,:,gg)   = tmptmpdata;
                    end;
                end;
                leg = opt.groupnames;
            else tmpdata = real(data{c,g}); 
                leg = {};;
            end;
            if ~isempty(opt.filter), tmpdata = myfilt(tmpdata, 1000/(allx(2)-allx(1)), 0, opt.filter); end;
            
            % plotting options
            % ----------------
            plotopt = { allx };
            if strcmpi(opt.plottopo, 'on'),
                if strcmpi(opt.plotsubjects, 'off') & strcmpi(opt.singlesubject, 'off') tmpdata = squeeze(real(nan_mean(tmpdata,3))); end;
            elseif strcmpi(opt.plotsubjects, 'off') & strcmpi(opt.singlesubject, 'off') tmpdata = squeeze(real(nan_mean(tmpdata,2))); 
            end;
            tmpdata = squeeze(permute(tmpdata, [2 1 3]));
            if strcmpi(opt.plottopo, 'on'), highlight = 'background'; else highlight = 'bottom'; end;
            if strcmpi(opt.plotgroups, 'together') &  isempty(opt.condstats) & ...
                             ~isnan(opt.threshold) & ~isempty(opt.groupstats)
                plotopt = { plotopt{:} 'maskarray' };
                tmpdata = { tmpdata pgroupplot{c}' };
            elseif strcmpi(opt.plotconditions, 'together') &  isempty(opt.groupstats) & ...
                                     ~isnan(opt.threshold) & ~isempty(opt.condstats)
                plotopt = { plotopt{:} 'maskarray' };
                tmpdata = { tmpdata pcondplot{c}' };
            end;
            plotopt = { plotopt{:} 'highlightmode', highlight };
            if strcmpi(opt.plotsubjects, 'on')
                plotopt = { plotopt{:} 'plotmean' 'on' 'plotindiv' 'on' };
            else
                plotopt = { plotopt{:} 'plotmean' 'off' };
            end;
            plotopt = { plotopt{:} 'ylim' opt.ylim 'legend' leg };
            plotopt = { plotopt{:}  'xlabel' xlab 'ylabel' ylab };
            
            % plot
            % ----
            if strcmpi(opt.figure, 'on'), tmpcol = col; else tmpcol = col(mod(c*ncplot+g,length(col))+1); end;
            if strcmpi(opt.plottopo, 'on'), 
                metaplottopo(tmpdata, 'chanlocs', opt.chanlocs, 'plotfunc', 'plotcurve', ...
                    'plotargs', { plotopt{:} }, 'datapos', [2 3]);
            elseif iscell(tmpdata)
                 plotcurve( allx, tmpdata{1}, 'colors', tmpcol, 'maskarray', tmpdata{2}, plotopt{3:end}); xlabel(xlab); ylabel(ylab);
            else plotcurve( allx, tmpdata, 'colors', tmpcol, plotopt{2:end});
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
        if g == ngplot & ng > 1 & ~isempty(opt.groupstats)            
            if ~strcmpi(opt.plotgroups, 'together') | ~isempty(opt.condstats) | isnan(opt.threshold)
                if strcmpi(opt.plotgroups, 'together'),         mysubplot(ncplot+addr, ngplot+addc, 2 + (c-1)*(ngplot+addc), opt.subplot); ci = g;
                elseif strcmpi(opt.plotconditions, 'together'), mysubplot(ncplot+addr, ngplot+addc, ngplot + 1, opt.subplot); ci = c;
                else                                            mysubplot(ncplot+addr, ngplot+addc, ngplot + 1 + (c-1)*(ngplot+addc), opt.subplot); ci = 1;
                end;
                if strcmpi(opt.plotconditions, 'together'), condnames = 'Conditions'; else condnames = opt.condnames{c}; end;
                if ~isnan(opt.threshold)
                     tmptitle = sprintf('%s (p<%.4f)', condnames, opt.threshold);
                     if strcmpi(opt.plottopo, 'on'), 
                          metaplottopo({zeros(size(pgroupplot{c}')) pgroupplot{c}'}, 'chanlocs', opt.chanlocs, 'plotfunc', 'plotcurve', ...
                              'plotargs', { allx 'maskarray' statopt{:} }, 'datapos', [2 3], 'title', tmptitle);
                     else plotcurve(allx, zeros(size(allx)), 'maskarray', mean(pgroupplot{c},2), 'ylim', [0.1 1], 'title', tmptitle, statopt{:});
                     end;
                else
                     if strcmpi(opt.plottopo, 'on'), 
                          metaplottopo(pgroupplot{c}', 'chanlocs', opt.chanlocs, 'plotfunc', 'plotcurve', ...
                              'plotargs', { allx statopt{:} }, 'datapos', [2 3], 'title', [ condnames ' (p-value)' ]);
                     else plotcurve(allx, mean(pgroupplot{c},2), 'title', [ condnames ' (p-value)' ], statopt{:});
                     end;
                end;
            end;
        end;
    end;
end;

for g = 1:ng
    % statistics accross conditions
    % -----------------------------
    if ~isempty(opt.condstats) & nc > 1
        if ~strcmpi(opt.plotconditions, 'together') | ~isempty(opt.groupstats) | isnan(opt.threshold)
            if strcmpi(opt.plotgroups, 'together'),         mysubplot(ncplot+addr, ngplot+addc, 1 + c*(ngplot+addc), opt.subplot); ci = g;
            elseif strcmpi(opt.plotconditions, 'together'), mysubplot(ncplot+addr, ngplot+addc, g + ngplot+addc, opt.subplot); ci = c;
            else                                            mysubplot(ncplot+addr, ngplot+addc, g + c*(ngplot+addc), opt.subplot); ci = 1;
            end;
            if strcmpi(opt.plotgroups, 'together'), groupnames = 'Groups'; else groupnames = opt.groupnames{g}; end;
            if ~isnan(opt.threshold)
                 tmptitle = sprintf('%s (p<%.4f)', groupnames, opt.threshold);
                 if strcmpi(opt.plottopo, 'on'), 
                      metaplottopo({zeros(size(pcondplot{g}')) pcondplot{g}'}, 'chanlocs', opt.chanlocs, 'plotfunc', 'plotcurve', ...
                          'plotargs', { allx 'maskarray' statopt{:} }, 'datapos', [2 3], 'title', tmptitle);
                 else plotcurve(allx, zeros(size(allx)), 'maskarray', mean(pcondplot{g},2), 'ylim', [0.1 1], 'title', tmptitle, statopt{:});
                 end;
            else
                 if strcmpi(opt.plottopo, 'on'), 
                      metaplottopo(pcondplot{g}', 'chanlocs', opt.chanlocs, 'plotfunc', 'plotcurve', ...
                          'plotargs', { allx statopt{:} }, 'datapos', [2 3], 'title', [ groupnames ' (p-value)' ]);
                 else plotcurve(allx, mean(pcondplot{g},2), 'title',  [ groupnames ' (p-value)' ],statopt{:});
                 end;
            end;
        end;
    end;
end;

% statistics accross group and conditions
% ---------------------------------------
if ~isempty(opt.groupstats) & ~isempty(opt.condstats) & ng > 1 & nc > 1
    mysubplot(ncplot+addr, ngplot+addc, ngplot + 1 + ncplot*(ngplot+addr), opt.subplot);
    if ~isnan(opt.threshold)
         tmptitle = sprintf('Interaction (p<%.4f)', opt.threshold);
         if strcmpi(opt.plottopo, 'on'), 
              metaplottopo({zeros(size(pinterplot')) pinterplot'}, 'chanlocs', opt.chanlocs, 'plotfunc', 'plotcurve', ...
                  'plotargs', { allx 'maskarray' statopt{:} }, 'datapos', [2 3], 'title', tmptitle);
         else plotcurve(allx, zeros(size(allx)), 'maskarray', mean(pinterplot,2), 'ylim', [0.1 1], 'title', tmptitle, statopt{:});
              xlabel(xlab); ylabel('-log10(p)');
        end;
    else
         if strcmpi(opt.plottopo, 'on'), 
              metaplottopo(pinterplot', 'chanlocs', opt.chanlocs, 'plotfunc', 'plotcurve', ...
                  'plotargs', { allx statopt{:} }, 'datapos', [2 3], 'title', 'Interaction (p-value)');
         else plotcurve(allx, mean(pinterplot,2), 'title', 'Interaction (p-value)', statopt{:});
         end;
    end;
end;  

% axis limit and legend
% ---------------------
for c = 1:ncplot
    for g = 1:ngplot
        if isempty(opt.ylim) & strcmpi(opt.plottopo, 'off')
            set(hdl(c,g), 'ylim', tmplim);
        end;
        if strcmpi(opt.plotgroups, 'together'),                        fig_title = opt.condnames{c};
        elseif strcmpi(opt.plotconditions, 'together'),                fig_title = opt.groupnames{g};
        elseif isempty(opt.condnames{c}) | isempty(opt.groupnames{g}), fig_title = [ opt.condnames{c} opt.groupnames{g} ];
        else                                                           fig_title = [ opt.condnames{c} ', ' opt.groupnames{g} ];
        end;
        if ~isempty(opt.compinds), fig_title = [ 'Comp.' int2str(opt.compinds{c,g}) ', ' fig_title ]; end;
        if ~isempty(opt.subject),  fig_title = [ opt.subject ', ' fig_title ]; end;
        set(get(hdl(c,g), 'title'), 'string', fig_title);
    end;
end;

if strcmpi(opt.plottopo, 'off'), 
    axcopy;
    % remove axis labels (for most but not all)
    % ------------------
    if strcmpi(opt.subplot, 'transpose')
        for c = 1:size(hdl,2)
            for g = 1:size(hdl,1)
                axes(hdl(g,c));
                if c ~= 1 & size(hdl,2) ~=1, xlabel(''); legend off; end;
                if g ~= 1 & size(hdl,1) ~= 1, ylabel(''); legend off; end;
            end;
        end;
    else
        for c = 1:size(hdl,1)
            for g = 1:size(hdl,2)
                axes(hdl(c,g));
                if g ~= 1 & size(hdl,2) ~=1, ylabel(''); legend off; end;
                if c ~= size(hdl,1) & size(hdl,1) ~= 1, xlabel(''); legend off; end;
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
