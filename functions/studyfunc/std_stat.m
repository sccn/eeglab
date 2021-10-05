% std_stat() - compute statistics for ERP/spectral traces or ERSP/ITC images
%              of a component or channel cluster in a STUDY.
% Usage:
%          >> [pcond, pgroup, pinter, statscond, statsgroup, statsinter] = std_stat( data, 'key', 'val', ...)
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
%                  abbreviations are still functional.
%  'naccu'       - [integer] Number of surrogate averages to accumulate when
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
%                  uses Fieldtrip statistical functions. Default is 'eeglab'.
%
% Fieldtrip statistics options:
%  'fieldtripnaccu'   - 'numrandomization' Fieldtrip parameter
%  'fieldtripalpha'   - 'alpha' Fieldtrip parameter. Default is 0.05.
%  'fieldtripmethod'  - 'method' Fieldtrip parameter. Default is 'analytic'
%  'fieldtripmcorrect' - 'mcorrect' Fieldtrip parameter. Default is 'none'.
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

function [pcond, pgroup, pinter, statscond, statsgroup, statsinter] = std_stat(data, varargin)

pgroup = {};
pcond  = {};
pinter = {};
if nargin < 1
    help std_stat;
    return;
end

% decode inputs
% -------------
if ~isempty(varargin) && isstruct(varargin{1})
    opt   = varargin{1};
    varargin(1) = [];
else
    opt = [];
end
if ~isempty(varargin) ||isempty(opt);
    opt = pop_statparams(opt, varargin{:});
end

if ~isfield(opt, 'paired'), opt.paired = { 'off' 'off' }; end
if ~isnan(opt.eeglab.alpha(1)) && isempty(opt.eeglab.naccu), opt.eeglab.naccu = 1/opt.eeglab.alpha(end)*2; end
if any(any(cellfun('size', data, 2)==1)), opt.groupstats = 'off'; opt.condstats = 'off'; end
if strcmpi(opt.eeglab.mcorrect, 'fdr'), opt.eeglab.naccu = opt.eeglab.naccu*20; end
if isempty(opt.eeglab.naccu), opt.eeglab.naccu = 2000; end
if ~isreal(data{1})
    fprintf('*** ITC significance - converting complex values to absolute amplitude ***\n');
    for ind = 1:length(data(:))
        data{ind} = abs(data{ind});
    end
end
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
            [F, df, pval] = statcond(data(:,g), 'method', opt.eeglab.method, 'naccu', opt.eeglab.naccu, 'paired', opt.paired{1});
            pcond{g}     = squeeze(pval);
            statscond{g} = squeeze(F);
        end
    end
    if strcmpi(opt.groupstats, 'on') && ng > 1
        for c = 1:nc
            [F, df, pval] = statcond(data(c,:), 'method', opt.eeglab.method, 'naccu', opt.eeglab.naccu, 'paired', opt.paired{2});
            pgroup{c}     = squeeze(pval);
            statsgroup{c} = squeeze(F);
        end
    end
    if ( strcmpi(opt.groupstats, 'on') || strcmpi(opt.condstats, 'on') ) && ng > 1 && nc > 1
        opt.paired = sort(opt.paired); % put 'off' first if present
        [F, df, pval] = statcond(data, 'method', opt.eeglab.method, 'naccu', opt.eeglab.naccu, 'paired', opt.paired{1});
        for index = 1:length(pval)
            pinter{index} = squeeze(pval{index});
            statsinter{index} = squeeze(F{index});
        end
    end
    
    if ~isempty(opt.groupstats) || ~isempty(opt.condstats)
        if ~strcmpi(opt.eeglab.mcorrect, 'none'),
            disp([ 'Applying ' upper(opt.eeglab.mcorrect) ' correction for multiple comparisons' ]);
            for ind = 1:length(pcond),  pcond{ind}  = mcorrect( pcond{ind} , opt.eeglab.mcorrect ); end
            for ind = 1:length(pgroup), pgroup{ind} = mcorrect( pgroup{ind}, opt.eeglab.mcorrect ); end
            if ~isempty(pinter),
                pinter{1} = mcorrect(pinter{1}, opt.eeglab.mcorrect);
                pinter{2} = mcorrect(pinter{2}, opt.eeglab.mcorrect);
                pinter{3} = mcorrect(pinter{3}, opt.eeglab.mcorrect);
            end
        end
        if ~isnan(opt.eeglab.alpha)
            for ind = 1:length(pcond),  pcond{ind}  = applythreshold(pcond{ind},  opt.eeglab.alpha); end
            for ind = 1:length(pgroup), pgroup{ind} = applythreshold(pgroup{ind}, opt.eeglab.alpha); end
            for ind = 1:length(pinter), pinter{ind} = applythreshold(pinter{ind}, opt.eeglab.alpha); end
        end
    end
else
    if ~exist('ft_freqstatistics'), error('Install Fieldtrip-lite to use Fieldtrip statistics'); end

    % Fieldtrip statistics
    % --------------------
    params = {};
    if strcmpi(opt.fieldtrip.mcorrect, 'cluster')
        params = eval( [ '{' opt.fieldtrip.clusterparam '}' ]);
        if isempty(opt.fieldtrip.channelneighbor), opt.fieldtrip.channelneighbor = struct([]); end
        params = { params{:} 'neighbours' opt.fieldtrip.channelneighbor }; % channelneighbor is empty if only one channel selected
    end
    params = { params{:} 'method', opt.fieldtrip.method, 'naccu', opt.fieldtrip.naccu 'mcorrect' opt.fieldtrip.mcorrect 'alpha' opt.fieldtrip.alpha 'numrandomization' opt.fieldtrip.naccu };
    params = { params{:} 'structoutput' 'on' }; % before if ~isnan(opt.fieldtrip.alpha), end
    
    pinter = cell(1,3);
    statsinter = cell(1,3);
    if strcmpi(opt.condstats, 'on') && nc > 1
        newdata = data(:,1);
        for g = 1:ng
            % marginal effect
            [F, df, pval]  = statcondfieldtrip(data(:,g), 'paired', opt.paired{1}, params{:});
            pcond{g}     = applymask( F, opt.fieldtrip);
            statscond{g} = squeeze(F.stat);
            
            % concatenate data for main statistics
            if g > 1 && ng > 1
                for c = 1:nc
                    switch ndims(data{1,1})
                        case 2, newdata{c,1}(:,end+1:end+size(data{c,g},2)) = data{c,g};
                        case 3, newdata{c,1}(:,:,end+1:end+size(data{c,g},3)) = data{c,g};
                        case 4, newdata{c,1}(:,:,:,end+1:end+size(data{c,g},4)) = data{c,g};
                    end
                end
            end
        end
        
        % main statistics
        if ng > 1
            [F, df, pval] = statcondfieldtrip(newdata, 'paired', opt.paired{1}, params{:});
            pinter{1}     = applymask(F, opt.fieldtrip);
            statsinter{1} = squeeze(F.stat);
        end
    else
        pcond = {};
    end
    if strcmpi(opt.groupstats, 'on') && ng > 1
        newdata = data(1,:);
        for c = 1:nc
            % marginal effect
            [F, df, pval]   = statcondfieldtrip(data(c,:), 'paired', opt.paired{2}, params{:});
            pgroup{c}     = applymask( F, opt.fieldtrip);
            statsgroup{c} = squeeze(F.stat);
            
            % concatenate data for main statistics
            if c > 1 && nc > 1
                for g = 1:ng
                    switch ndims(data{1,1})
                        case 2, newdata{1,g}(:,end+1:end+size(data{c,g},2)) = data{c,g};
                        case 3, newdata{1,g}(:,:,end+1:end+size(data{c,g},3)) = data{c,g};
                        case 4, newdata{1,g}(:,:,:,end+1:end+size(data{c,g},4)) = data{c,g};
                    end
                end
            end
        end
        
        % main statistics
        if nc > 1
            [F, df, pval] = statcondfieldtrip(newdata, 'paired', opt.paired{1}, params{:});
            pinter{2}     = applymask(F, opt.fieldtrip);
            statsinter{2} = squeeze(F.stat);
        end
    else
        pgroup = {};
    end
end

% apply mask for fieldtrip data
% -----------------------------
function p = applymask(F, fieldtrip)

if ~isnan(fieldtrip.alpha), p = squeeze(F.mask);
else
    p = squeeze(F.pval);
    if ~strcmpi(fieldtrip.mcorrect, 'none')
        p(~F.mask) = 1;
    end
end

% apply stat threshold to data for EEGLAB stats
% ---------------------------------------------
function newdata = applythreshold(data, threshold)

threshold = sort(threshold);
newdata = zeros(size(data));
for index = 1:length(threshold)
    inds = data < threshold(index);
    data(inds)    = 1;
    newdata(inds) = length(threshold)-index+1;
end

% compute correction for multiple comparisons
% -------------------------------------------
function pvals = mcorrect(pvals, method);

switch method
    case {'no' 'none'}, return;
    case 'bonferoni', pvals = pvals*prod(size(pvals));
    case 'holms',     [tmp ind] = sort(pvals(:)); [tmp ind2] = sort(ind); pvals(:) = pvals(:).*(prod(size(pvals))-ind2+1);
    case 'fdr',       pvals = fdr(pvals);
    otherwise error(['Unknown method ''' method ''' for correction for multiple comparisons' ]);
end        
        
