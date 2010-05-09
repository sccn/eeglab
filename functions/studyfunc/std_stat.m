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
%  'threshold'   - [NaN|p-value] threshold for computing p-value. In this 
%                  function, it is only used to compute naccu above. NaN
%                  means that no threshold has been set. 
%  'mcorrect'    - ['none'|'fdr'] apply correcting for multiple comparisons.
%
% Outputs:
%  pcond        - [cell] condition pvalues or mask (0 or 1). One element per group.
%  pgroup       - [cell] group pvalues or mask (0 or 1). One element per condition.
%  pinter       - [cell] three elements, condition pvalues (group pooled),
%                 group pvalues (condition pooled) and interaction pvalues.
%
% Author: Arnaud Delorme, CERCO, CNRS, 2006-
% 
% See also: statcond()

%  'statmode'    - ['subjects'|'trials'] standard statistics are 
%                  'subjects' where the statistics is performed accross
%                  the mean ERSP (or ITC) of single subjects. For 'trials'
%                  statistics, the single-trial data epochs of all subjects
%                  are pooled together. This requires that they were
%                  saved on disk using option 'savetrials', 'on' at the time
%                  of computation. Note that these single-trial data
%                  may use several GB of disk space and that computation 
%                  of 'trials' statistics requires a lot of RAM.

% $Log: std_stat.m,v $
% Revision 1.12  2010/03/09 06:14:16  arno
% fixed default output
%
% Revision 1.11  2010/03/05 01:25:19  arno
% Fix threshold and interstat plotting
%
% Revision 1.10  2010/02/09 06:07:27  arno
% Fixed new title problem and implemented 3-level significance
%
% Revision 1.9  2010/01/28 20:46:01  arno
% fixed computing FDR for Anova interactions
%
% Revision 1.8  2009/08/29 04:24:56  arno
% new statistics
%
% Revision 1.6  2009/08/29 00:38:32  arno
% move all statistics to std_stat
%
% Revision 1.5  2009/08/04 23:21:25  arno
% naccu
%
% Revision 1.4  2009/05/31 02:22:10  arno
% Adding FDR and bootstrap to all STUDY functions
%
% Revision 1.3  2007/09/11 10:55:48  arno
% disable stats if not enough info
%
% Revision 1.2  2007/04/30 19:52:22  arno
% thrshold entry\
%
% Revision 1.1  2007/01/26 18:09:02  arno
% Initial revision
%

function [pcond, pgroup, pinter] = std_stat(data, varargin)

pgroup = {};
pcond  = {};
pinter = {};
if nargin < 1
    help std_stat;
    return;
end;

% decode inputs
% -------------
opt = {};
if nargin > 1
    if isstruct(varargin{1})
        g = struct(varargin{2:end});
        if isfield(g, 'condstats'),  varargin{1} = rmfield(varargin{1}, 'condstats'); end;
        if isfield(g, 'groupstats'), varargin{1} = rmfield(varargin{1}, 'groupstats'); end;
        tmpparams      = fieldnames(varargin{1}); tmpparams = tmpparams';
        tmpparams(2,:) = struct2cell(varargin{1});
        varargin       = { tmpparams{:} varargin{2:end} };
    end;
end;
opt = finputcheck( varargin, { 'threshold'   'real'    []               NaN;
                               'mcorrect'    'string'  { 'none' 'fdr' } 'none';
                               'naccu'       'integer' []               [];
                               'groupstats'  'string'  { 'on' 'off' }   'off';
                               'paired'      'cell'    { 'on' 'off' }   { 'on' 'on' };
                               'condstats'   'string'  { 'on' 'off' }   'off';
                               'statistics'  'string'  { 'param' 'perm' 'bootstrap' }       'param' }, ...
                               'std_stat', 'ignore');

if isstr(opt), error(opt); end;
if ~isnan(opt.threshold(1)) && isempty(opt.naccu), opt.naccu = 1/opt.threshold(end)*2; end;
if any(any(cellfun('size', data, 2)==1)), opt.groupstats = 'off'; opt.condstats = 'off'; end;
if strcmpi(opt.mcorrect, 'fdr'), opt.naccu = opt.naccu*20; end;
if isempty(opt.naccu), opt.naccu = 2000; end;
if strcmpi(opt.statistics, 'param') && ~isreal(data{1})
    fprintf('*** Cannot use parametric statistics for single-trial ITC significance ***\n');
    return;
end;
nc = size(data,1);
ng = size(data,2);

% compute significance mask
% --------------------------
if strcmpi(opt.condstats, 'on') && nc > 1
    for g = 1:ng
        [F df pval] = statcond(data(:,g), 'mode', opt.statistics, 'naccu', opt.naccu, 'paired', opt.paired{1}); 
        pcond{g} = squeeze(pval);
    end;
else
    pcond = {};
end;
if strcmpi(opt.groupstats, 'on') && ng > 1
    for c = 1:nc
        [F df pval] = statcond(data(c,:), 'mode', opt.statistics, 'naccu', opt.naccu, 'paired', opt.paired{2}); 
        pgroup{c} = squeeze(pval);
    end;
else
    pgroup = {};
end;
if ( strcmpi(opt.groupstats, 'on') && strcmpi(opt.condstats, 'on') ) & ng > 1 & nc > 1
    opt.paired = sort(opt.paired); % put 'off' first if present
    [F df pval] = statcond(data, 'mode', opt.statistics, 'naccu', opt.naccu, 'paired', opt.paired{1});
    for index = 1:length(pval)
        pinter{index} = squeeze(pval{index});
    end;
else
    pinter = {};
end;

if ~isempty(opt.groupstats) || ~isempty(opt.condstats)   
    if strcmpi(opt.mcorrect, 'fdr'), 
        disp('Applying FDR correction for multiple comparisons');
        for ind = 1:length(pcond),  pcond{ind} = fdr( pcond{ind} ); end;
        for ind = 1:length(pgroup), pgroup{ind} = fdr( pgroup{ind} ); end;
        if ~isempty(pinter), 
            pinter{1} = fdr(pinter{1}); 
            pinter{2} = fdr(pinter{2}); 
            pinter{3} = fdr(pinter{3}); 
        end;
    end;
    if ~isnan(opt.threshold)
        for ind = 1:length(pcond),  pcond{ind}  = applythreshold(pcond{ind},  opt.threshold); end;
        for ind = 1:length(pgroup), pgroup{ind} = applythreshold(pgroup{ind}, opt.threshold); end;
        for ind = 1:length(pinter), pinter{ind} = applythreshold(pinter{ind}, opt.threshold); end;
    end;
end;

function newdata = applythreshold(data, threshold)
    threshold = sort(threshold);
    newdata = zeros(size(data));
    for index = 1:length(threshold)
        inds = data < threshold(index);
        data(inds)    = 1;
        newdata(inds) = length(threshold)-index+1;
    end;
    