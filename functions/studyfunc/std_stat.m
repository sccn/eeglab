% std_stat() - compute statistics for ERP/spectral traces or ERSP/ITC images
%              of a component or channel cluster in a STUDY.
% Usage:
%          >> std_stat( data, 'key', 'val', ...)
% Inputs:
%  data  -  [cell array] mean data for each subject group and/or data
%           condition. For example, to compute mean ERPs statistics from a
%           STUDY for epochs of 800 frames in two conditions from three
%           groups of 12 subjects:
%
%           >> data = { [800x12] [800x12] [800x12];... % 3 groups, cond 1
%                       [800x12] [800x12] [800x12] };  % 3 groups, cond 2
%           >> pcond = std_stat(data, 'condstats', 'on');
%
%           By default, parametric statistics are computed across subjects
%           in the three groups. See below and >> help statcond
%           for more information about the statistical computations.
%
% Statistics options (EEGLAB):
%  'groupstats'  - ['on'|'off'] Compute (or not) statistics across groups.
%                  {default: 'off'}
%  'condstats'   - ['on'|'off'] Compute (or not) statistics across groups.
%                  {default: 'off'}
%  'method'      - ['parametric'|'permutation'] Type of statistics to use
%                  default is 'parametric'. 'perm' and 'param' legacy
%                  abreviations are still functional.
%  'naccu'       - [integer] Number of surrogate averages fo accumulate when
%                  computing permutation-based statistics. For example, to
%                  test p<0.01 use naccu>=200; for p<0.001, use naccu>=2000.
%                  If a non-NaN 'threshold' is set (see below) and 'naccu'
%                  is too low, it will be automatically increased. This
%                  keyword available only from the command line {default:500}
%  'alpha'       - [NaN|p-value] threshold for computing p-value. In this
%                  function, it is only used to compute naccu above. NaN
%                  means that no threshold has been set.
%  'mcorrect'    - ['none'|'fdr'] apply correcting for multiple comparisons.
%  'mode'        - ['eeglab'|'fieldtrip'] statistical framework to use. 
%                  'eeglab' uses EEGLAB statistical functions and 'fieldtrip'
%                  uses Fieldtrip statistical funcitons. Default is 'eeglab'.
%
% Fieldtrip statistics options:
%  'fieldtripnaccu'   - 'numrandomization' Fieldtrip parameter
%  'fieldtripalpha'   - 'alpha' Fieldtrip parameter. Default is 0.05.
%  'fieldtripmethod'  - 'method' Fieldtrip parameter. Default is 'analytic'
%  'fieldtripmcorrect' - 'mcorrect' Fieldtrip parameter. Default is 'no'.
%  'fieldtripclusterparam' - string or cell array for optional parameters
%                            for cluster correction method, see function
%                            ft_statistics_montecarlo for more information.
%  'fieldtripchannelneighbor' - Fieldtrip channel neighbour structure for 
%                               cluster correction method, see function
%                               std_prepare_neighbors for more information.

% Legacy parameters:
%   'threshold'  - now 'alpha'
%   'statistics' - now 'method'  
%
% Outputs:
%  pcond        - [cell] condition pvalues or mask (0 or 1) if an alpha value
%                 is selected. One element per group.
%  pgroup       - [cell] group pvalues or mask (0 or 1). One element per 
%                 condition.
%  pinter       - [cell] three elements, condition pvalues (group pooled),
%                 group pvalues (condition pooled) and interaction pvalues.
%  statcond     - [cell] condition statistic values (F or T).
%  statgroup    - [cell] group pvalues or mask (0 or 1). One element per 
%                 condition.
%  statinter    - [cell] three elements, condition statistics (group pooled),
%                 group statistics (condition pooled) and interaction F statistics.
%
% Author: Arnaud Delorme, CERCO, CNRS, 2006-
%
% See also: statcond()

% Copyright (C) Arnaud Delorme
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

function [pcond, pgroup, pinter, statscond, statsgroup, statsinter] = std_stat(data, varargin)

pgroup = {};
pcond  = {};
pinter = {};
if nargin < 1
    help std_stat;
    return;
end;

% decode inputs
% -------------
if isstruct(varargin{1})
    opt   = varargin{1};
else
    statstruct.etc = STUDY.etc;
    statstruct     = pop_statparams(statstruct, varargin{:});
    opt            = statstruct.etc.statistics;
end;

if ~isnan(opt.eeglab.alpha(1)) && isempty(opt.eeglab.naccu), opt.eeglab.naccu = 1/opt.eeglab.alpha(end)*2; end;
if any(any(cellfun('size', data, 2)==1)), opt.groupstats = 'off'; opt.condstats = 'off'; end;
if strcmpi(opt.eeglab.mcorrect, 'fdr'), opt.eeglab.naccu = opt.eeglab.naccu*20; end;
if isempty(opt.eeglab.naccu), opt.eeglab.naccu = 2000; end;
if strcmpi(opt.eeglab.method(1:3), 'par') && ~isreal(data{1})
    fprintf('*** Cannot use parametric statistics for single-trial ITC significance ***\n');
    return;
end;
nc = size(data,1);
ng = size(data,2);

% compute significance mask
% -------------------------
pcond  = {};
pgroup = {};
pinter = {};
statscond  = {};
statsgroup = {};
statsinter = {};
if strcmpi(opt.mode, 'eeglab')
    % EEGLAB statistics
    % -----------------
    if strcmpi(opt.condstats, 'on') && nc > 1
        for g = 1:ng
            [F df pval] = statcond(data(:,g), 'method', opt.eeglab.method, 'naccu', opt.eeglab.naccu, 'paired', opt.paired{1});
            pcond{g}     = squeeze(pval);
            statscond{g} = squeeze(F);
        end;
    end;
    if strcmpi(opt.groupstats, 'on') && ng > 1
        for c = 1:nc
            [F df pval] = statcond(data(c,:), 'method', opt.eeglab.method, 'naccu', opt.eeglab.naccu, 'paired', opt.paired{2});
            pgroup{c}     = squeeze(pval);
            statsgroup{c} = squeeze(F);
        end;
    else
    end;
    if ( strcmpi(opt.groupstats, 'on') && strcmpi(opt.condstats, 'on') ) & ng > 1 & nc > 1
        opt.paired = sort(opt.paired); % put 'off' first if present
        [F df pval] = statcond(data, 'method', opt.eeglab.method, 'naccu', opt.eeglab.naccu, 'paired', opt.paired{1});
        for index = 1:length(pval)
            pinter{index} = squeeze(pval{index});
            statsinter{index} = squeeze(F{index});
        end;
    end;
    
    if ~isempty(opt.groupstats) || ~isempty(opt.condstats)
        if strcmpi(opt.eeglab.mcorrect, 'fdr'),
            disp('Applying FDR correction for multiple comparisons');
            for ind = 1:length(pcond),  pcond{ind}  = mcorrect( pcond{ind} , opt.eeglab.mcorrect ); end;
            for ind = 1:length(pgroup), pgroup{ind} = mcorrect( pgroup{ind}, opt.eeglab.mcorrect ); end;
            if ~isempty(pinter),
                pinter{1} = mcorrect(pinter{1}, opt.eeglab.mcorrect);
                pinter{2} = mcorrect(pinter{2}, opt.eeglab.mcorrect);
                pinter{3} = mcorrect(pinter{3}, opt.eeglab.mcorrect);
            end;
        end;
        if ~isnan(opt.eeglab.alpha)
            for ind = 1:length(pcond),  pcond{ind}  = applythreshold(pcond{ind},  opt.eeglab.alpha); end;
            for ind = 1:length(pgroup), pgroup{ind} = applythreshold(pgroup{ind}, opt.eeglab.alpha); end;
            for ind = 1:length(pinter), pinter{ind} = applythreshold(pinter{ind}, opt.eeglab.alpha); end;
        end;
    end;
else
    % Fieldtrip statistics
    % --------------------
    params = {};
    if strcmpi(opt.fieldtrip.mcorrect, 'cluster')
        params = eval( [ '{' opt.fieldtrip.clusterparam '}' ]);
        if isempty(opt.fieldtrip.channelneighbor), opt.fieldtrip.channelneighbor = struct([]); end;
        params = { params{:} 'neighbours' opt.fieldtrip.channelneighbor };
    end;
    params = { params{:} 'method', opt.fieldtrip.method, 'naccu', opt.fieldtrip.naccu 'mcorrect' opt.fieldtrip.mcorrect 'alpha' opt.fieldtrip.alpha 'numrandomization' opt.fieldtrip.naccu };
        
    if strcmpi(opt.condstats, 'on') && nc > 1
        for g = 1:ng
            [F df pval] = statcondfieldtrip(data(:,g), 'paired', opt.paired{1}, params{:});
            pcond{g} = squeeze(pval);
            statscond{g} = squeeze(F);
        end;
    else
        pcond = {};
    end;
    if strcmpi(opt.groupstats, 'on') && ng > 1
        for c = 1:nc
            [F df pval] = statcondfieldtrip(data(c,:), 'paired', opt.paired{2}, params{:});
            pgroup{c} = squeeze(pval);
            statsgroup{c} = squeeze(F);
        end;
    else
        pgroup = {};
    end;
    if ( strcmpi(opt.groupstats, 'on') && strcmpi(opt.condstats, 'on') ) & ng > 1 & nc > 1
        opt.paired = sort(opt.paired); % put 'off' first if present
        [F df pval] = statcondfieldtrip(data, 'paired', opt.paired{1}, params{:});
        for index = 1:length(pval)
            pinter{index} = squeeze(pval{index});
            statsinter{index} = squeeze(F{index});
        end;
    else
        pinter = {};
    end;
end;

% apply stat threshold to data
% ----------------------------
function newdata = applythreshold(data, threshold)

threshold = sort(threshold);
newdata = zeros(size(data));
for index = 1:length(threshold)
    inds = data < threshold(index);
    data(inds)    = 1;
    newdata(inds) = length(threshold)-index+1;
end;

% compute correction for multiple comparisons
% -------------------------------------------
function pvals = mcorrect(pvals, method);

switch method
    case 'none', return;
    case 'bonferoni', pvals = pvals*prod(size(pvals));
    case 'holms',     [tmp ind] = sort(pvals(:)); [tmp ind2] = sort(ind); pvals(:) = pvals(:).*(prod(size(pvals))-ind2+1);
    case 'fdr',       pvals = fdr(pvals);
    otherwise error(['Unknown method ''' method ''' for correction for multiple comparisons' ]);
end;        
        
