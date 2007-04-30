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
%
% Outputs:
%  pcond        - [cell] condition pvalues. One element per group.
%  pgroup       - [cell] group pvalues. One element per condition.
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

% $Log: not supported by cvs2svn $
% Revision 1.1  2007/01/26 18:09:02  arno
% Initial revision
%

function [pcond, pgroup, pinter] = std_stat(data, varargin)

pgroup = [];
pcond  = [];
pinter = [];
if nargin < 1
    help std_plot;
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
opt = finputcheck( varargin, { 'threshold'   'real'   []              NaN;
                               'naccu'       'integer' []             [];
                               'groupstats'   'string' { 'on' 'off' }   'off';
                               'condstats'    'string' { 'on' 'off' }   'off';
                               'statistics'  'string' { 'param' 'perm' }       'param' }, ...
                               'std_stat', 'ignore');
if isstr(opt), error(opt); end;
if ~isnan(opt.threshold) & isempty(opt.naccu), opt.naccu = 1/opt.threshold*2; end;
if isempty(opt.naccu), opt.naccu = 500; end;
nc = size(data,1);
ng = size(data,2);

% compute significance mask
% --------------------------
if strcmpi(opt.condstats, 'on') & nc > 1
    for g = 1:ng
        [F df pval] = statcond(data(:,g), 'mode', opt.statistics, 'naccu', opt.naccu); 
        pcond{g} = squeeze(pval);
    end;
else
    pcond = {};
end;
if strcmpi(opt.groupstats, 'on') & ng > 1
    for c = 1:nc
        [F df pval] = statcond(data(c,:), 'mode', opt.statistics, 'naccu', opt.naccu); 
        pgroup{c} = squeeze(pval);
    end;
else
    pgroup = {};
end;
if ( strcmpi(opt.groupstats, 'on') & strcmpi(opt.condstats, 'on') ) & ng > 1 & nc > 1
    [F df pval] = statcond(data, 'mode', opt.statistics, 'naccu', opt.naccu);
    for index = 1:length(pval)
        pinter{index} = squeeze(pval{index});
    end;
else
    pinter = {};
end;
