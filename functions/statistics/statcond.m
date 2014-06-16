% statcond()  - compare two or more data conditions statistically using 
%               standard parametric or nonparametric permutation-based ANOVA 
%               (1-way or 2-way) or t-test methods. Parametric testing uses 
%               fcdf() from the Matlab Statistical Toolbox.
% Usage:
%          >> [stats, df, pvals, surrog] = statcond( data, 'key','val'... );
%
% Inputs:
%   data       = one-or two-dimensional cell array of data matrices. 
%                   For nonparametric, permutation-based testing, the 
%                last dimension of the data arrays (which may be of up to 
%                4 dimensions) is permuted across conditions, either in
%                a 'paired' fashion (not changing the, e.g., subject or
%                trial order in the last dimension) or in an umpaired 
%                fashion (not respecting this order). If the number of 
%                elements in the last dimension is not the same across 
%                conditions, the 'paired' option is turned 'off'.  Note: 
%                All other dimensions MUST be constant across conditions. 
%                   For example, consider a (1,3) cell array of matrices 
%                of size (100,20,x) each holding a (100,20) time/frequency 
%                transform from each of x subjects. Only the last dimension 
%                (here x, the number of subjects) may differ across the 
%                three conditions. 
%                   The test used depends on the size of the data array input.
%                When the data cell array has 2 columns and the data are 
%                paired, a paired t-test is performed; when the data are 
%                unpaired, an unpaired t-test is performed. If 'data' 
%                has only one row (paired or unpaired) and more than 2 
%                columns, a one-way ANOVA is performed. If the data cell 
%                array contains several rows and columns, and the data is
%                paired, a two-way repeated measure ANOVA is performed. 
%                NOTE THAT IF THE DATA is unpaired, EEGLAB will use a 
%                balanced 1 or 2 way ANOVA and parametric results might not 
%                be meaningful (bootstrap and permstatcondutation should be fine).
%
% Optional inputs:
%   'paired'   = ['on'|'off'] pair the data array {default: 'on' unless 
%                the last dimension of data array is of different lengths}.
%                For two independent variables, this input is a cell array,
%                for example { 'on' 'off' } indicating that the first
%                independent variable is paired and the second is not.
%   'method'   = ['perm'|'bootstrap'|'param'] method for computing the p-values:
%                 'param' or 'parametric' = parametric testing (standard ANOVA
%                                           or t-test); 
%                 'perm' or 'permutation' = non-parametric testing using 
%                                           surrogate data
%                 'bootstrap' = non-parametric bootstrap 
%                  made by permuting the input data {default: 'param'}
%   'naccu'    = [integer] Number of surrogate data copies to use in 'perm' 
%                 or 'bootstrap' method estimation (see above) {default: 200}.
%   'verbose'  = ['on'|'off'] print info on the command line {default: 'on'}.
%   'variance' = ['homegenous'|'inhomogenous'] this option is exclusively
%                for parametric statistics using unpaired t-test. It allows
%                to compute a more accurate value for the degree of freedom
%                using the formula for inhomogenous variance (see
%                ttest2_cell function). Default is 'inhomegenous'.
%   'surrog'   = surrogate data array (see output).
%   'stats'    = F- or T-value array (see output).
%   'tail'     = ['one'|'two'] run one-tailed (F-test) or two tailed
%                (T-test). This option is only relevant when using the
%                'surrog' input. Otherwise it is ignored.
%   'forceanova' = ['on'|'off'] force the use of ANOVA calculation even
%                for 2x1 designs. Default is 'off'.
%   'alpha'    = [float] p-value threshold value. Allow returning
%                confidence intervals and mask (requires structoutput below).
%   'structoutput' = ['on'|'off'] return an output structure instead of 
%                the regular output. Allow to output mask and confidence
%                intervals.
%
% Legacy parameters:
%   'threshold' - now 'alpha'
%   'mode'      - now 'method'
%
% Outputs:
%   stats      = F- or T-value array of the same size as input data without 
%                the last dimension. A T value is returned only when the data 
%                includes exactly two conditions.
%   df         = degrees of freedom, a (2,1) vector, when F-values are returned
%   pvals      = array of p-values. Same size as input data without the last
%                data dimension. All returned p-values are two-tailed.
%   surrog     = surrogate data array (same size as input data with the last 
%                dimension filled with a number ('naccu') of surrogate data sets.
%
% Important note: When a two-way ANOVA is performed, outputs are cell arrays
%                 with three elements: output(1) = row effects; 
%                 output(2) = column effects; output(3) = interactions
%                 between rows and columns.
%
% Examples:
%      >> a = { rand(1,10) rand(1,10)+0.5 }; % pseudo 'paired' data vectors
%         [t df pvals] = statcond(a);        % perform paired t-test
%           pvals =                  
%              5.2807e-04 % standard t-test probability value
%         % Note: for different rand() outputs, results will differ.
%
%         [t df pvals surog] = statcond(a, 'method', 'perm', 'naccu', 2000); 
%           pvals =
%              0.0065 % nonparametric t-test using 2000 permuted data sets
%
%         a = { rand(2,11) rand(2,10) rand(2,12)+0.5 }; % pseudo 'unpaired' 
%         [F df pvals] = statcond(a); % perform an unpaired ANOVA 
%           pvals =
%              0.00025 % p-values for difference between columns 
%              0.00002 % for each data row
%
%         a = { rand(3,4,10) rand(3,4,10) rand(3,4,10); ...
%               rand(3,4,10) rand(3,4,10) rand(3,4,10)+0.5 }; 
%         % pseudo (2,3)-condition data array, each entry containing 
%         %                                    ten (3,4) data matrices
%         [F df pvals] = statcond(a);  % perform a paired 2-way ANOVA 
%         % Output:
%           pvals{1} % a (3,4) matrix of p-values; effects across rows
%           pvals{2} % a (3,4) matrix of p-values; effects across colums 
%           pvals{3} % a (3,4) matrix of p-values; interaction effects
%                                      % across rows and columns
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 2005-
%         With thanks to Robert Oostenveld for fruitful discussions 
%         and advice on this function.
%
% See also: anova1_cell(), anova2_cell(), anova2rm_cell, fcdf()

% perform a paired t-test
% -----------------------
% a = { rand(2,10) rand(2,10) };
% [t df pval] = statcond(a); pval
% [h p t stat] = ttest( a{1}(1,:), a{2}(1,:)); p
% [h p t stat] = ttest( a{1}(2,:), a{2}(2,:)); p
%
% compare significance levels
% --------------------------
% a = { rand(1,10) rand(1,10) }; 
% [F df pval] = statcond(a, 'method', 'perm', 'naccu', 200); pval
% [h p t stat] = ttest( a{1}(1,:), a{2}(1,:)); p

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

function [ ori_vals, df, pvals, surrogval ] = statcond( data, varargin );
    
    if nargin < 1
        help statcond;
        return;
    end;
    try, warning('off', 'MATLAB:divideByZero'); catch, end;    
    
    if exist('finputcheck')
        g = finputcheck( varargin, { 'naccu'      'integer'   [1 Inf]             200;
                                     'method'     'string'    { 'param','parametric','perm','permutation','bootstrap' }  'param';
                                     'mode'       'string'    { }                 '';
                                     'paired'     'string'    { 'on','off' }      'on'; 
                                     'surrog'     { 'real','cell' }      []       []; 
                                     'stats'      { 'real','cell' }      []       []; 
                                     'structoutput' 'string'  { 'on','off' }      'off'; 
                                     'forceanova'   'string'  { 'on','off' }      'off'; 
                                     'arraycomp'  'string'    { 'on','off' }      'on'; 
                                     'alpha'      'real'      []                  NaN;
                                     'tail'       'string'    { 'one','both','upper','lower'}    'both'; 
                                     'variance'   'string'    { 'homogenous','inhomogenous' }    'inhomogenous'; 
                                     'returnresamplingarray' 'string'    { 'on','off' }      'off'; 
                                     'verbose'    'string'    { 'on','off' }      'on' }, 'statcond');
        if isstr(g), error(g); end;
    else
        g = struct(varargin{:});
        if ~isfield(g, 'naccu'),     g.naccu = 200; end;
        if ~isfield(g, 'method'),    g.method  = 'param'; end;
        if ~isfield(g, 'paired'),    g.paired = 'on'; end;
        if ~isfield(g, 'surrog'),    g.surrog = []; end;
        if ~isfield(g, 'orivals'),   g.orivals = []; end;
        if ~isfield(g, 'arraycomp'), g.arraycomp = 'on'; end;
        if ~isfield(g, 'verbose'),   g.verbose = 'on'; end;
        if ~isfield(g, 'tail'),      g.tail = 'both'; end;
        if ~isfield(g, 'variance'),  g.variance = 'homogenous'; end;
        if ~isfield(g, 'structoutput'), g.structoutput = 'on'; end;
        if ~isfield(g, 'returnresamplingarray'),   g.returnresamplingarray = 'off'; end;
    end;
    if ~isempty(g.mode), g.method = g.mode; end;
    
    if strcmpi(g.method, 'parametric'), g.method = 'param'; end;
    if strcmpi(g.method, 'permutation'), g.method = 'perm'; end;
    if strcmpi(g.verbose, 'on'), verb = 1; else verb = 0; end;
    if strcmp(g.method, 'param' ) && exist('fcdf') ~= 2
      myfprintf('on',['statcond(): parametric testing requires fcdf() \n' ...
               '            from the Matlab StatsticaL Toolbox.\n' ...
               '            Running nonparametric permutation tests\n.']);
      g.method = 'perm';
    end
    if size(data,2) == 1, data  = transpose(data); end; % cell array transpose
    g.naccu = round(g.naccu);
    
    % reshape matrices
    % ----------------
    nd = size(data{1});
    nd = nd(1:end-1);
    for index = 1:prod(size(data))
        data{index} = reshape(data{index}, [prod(nd) size(data{index},myndims(data{index}))]);
    end;    
    
    if ~strcmpi(g.method, 'param') && isempty(g.surrog)
         tmpsize   = size(data{1});
         surrogval = zeros([ tmpsize(1:end-1) g.naccu ], 'single');
    else surrogval = [];
    end;
    
    % check for NaNs or Inf
    % ---------------------
    for iDat = 1:length(data(:))
        if any(isnan(reshape(data{iDat}, prod(size(data{iDat})),1))) || ...
                any(isinf(reshape(data{iDat}, prod(size(data{iDat})),1)))
            error('Statcond: One of the input array contains NaNs or Infinite values');
        end;
    end;
        
    % bootstrap flag
    % --------------
    if strcmpi(g.method, 'bootstrap'), bootflag = 1;
    else                               bootflag = 0;
    end;
    
    if isempty(g.surrog)
        % test if data can be paired
        % --------------------------
        if length(unique(cellfun('size', data, ndims(data{1}) ))) > 1
            g.paired = 'off'; 
        end;
        if strcmpi(g.paired, 'on')
             pairflag = 1;
        else pairflag = 0;
        end;

        % return resampling array
        % -----------------------
        if strcmpi(g.returnresamplingarray, 'on')
            [ datavals datalen datadims ] = concatdata( data );
            if strcmpi(g.arraycomp, 'on')
                ori_vals = surrogdistrib( data, 'method', g.method, 'pairing', g.paired, 'naccu', g.naccu);
            else
                ori_vals = surrogdistrib( data, 'method', g.method, 'pairing', g.paired);
            end;
            return;
        end;
        
        % text output
        % -----------
        myfprintf(verb,'%d x %d, ', size(data,1), size(data,2));
        if strcmpi(g.paired, 'on')
             myfprintf(verb,'paired data, ');
        else myfprintf(verb,'unpaired data, ');
        end;
        if size(data,1) == 1 && size(data,2) == 2
             myfprintf(verb,'computing T values\n');
        else myfprintf(verb,'computing F values\n');
        end;
        if size(data,1) > 1 
            if strcmpi(g.paired, 'on')
                 myfprintf(verb,'Using 2-way repeated measure ANOVA\n');
            else myfprintf(verb,'Using balanced 2-way ANOVA (not suitable for parametric testing, only bootstrap)\n');
            end;
        elseif size(data,2) > 2
            if strcmpi(g.paired, 'on')
                 myfprintf(verb,'Using 1-way repeated measure ANOVA\n');
            else myfprintf(verb,'Using balanced 1-way ANOVA (equivalent to Matlab anova1)\n');
            end;
        else
            if strcmpi(g.paired, 'on')
                 myfprintf(verb,'Using paired t-test\n');
            else myfprintf(verb,'Using unpaired t-test\n');
            end;
        end;
        if ~strcmpi(g.method, 'param')
            if bootflag, myfprintf(verb,'Bootstraps (of %d):', g.naccu);
            else         myfprintf(verb,'Permutations (of %d):', g.naccu);
            end;
        end;
    end;
    
    tail = g.tail;
    if isempty(g.surrog)
        if size(data,1) == 1, % only one row

            if size(data,2) == 2 && strcmpi(g.forceanova, 'off')

                % paired t-test (very fast)
                % -------------
                [ori_vals df] = ttest_cell_select(data, g.paired, g.variance);

                if strcmpi(g.method, 'param')
                    
                    % Check if exist tcd.m file from the Statistics Toolbox (Bug 1352 )
                    if exist('tcdf','file') == 2  & license('test', 'Statistics_Toolbox')
                        pvals = 2*tcdf(-abs(ori_vals), df);
                    else
                        pvals = 2*mytcdf(-abs(ori_vals), df);
                    end
                    
                    pvals = reshape(pvals, size(pvals));
                else
                    if strcmpi(g.arraycomp, 'on')
                        try
                            myfprintf(verb,'...');
                            res = surrogdistrib( data, 'method', g.method, 'pairing', g.paired, 'naccu', g.naccu);
                            surrogval = ttest_cell_select( res, g.paired, g.variance);
                        catch,
                           lasterr
                           myfprintf(verb,'\nSuperfast array computation failed because of memory limitation, reverting to standard computation');
                           g.arraycomp = 'off';
                        end;
                    end;
                    if strcmpi(g.arraycomp, 'off')
                        [res precomp] = surrogdistrib( data, 'method', g.method, 'pairing', g.paired);
                        for index = 1:g.naccu
                            res = surrogdistrib( {}, 'precomp', precomp);
                            if mod(index, 10) == 0, myfprintf(verb,'%d ', index); end;
                            if mod(index, 100) == 0, myfprintf(verb,'\n'); end;
                            if myndims(res{1}) == 1
                                 surrogval(index)     = ttest_cell_select(res, g.paired, g.variance);
                            else surrogval(:,index)   = ttest_cell_select(res, g.paired, g.variance);
                            end;
                        end;
                    end;
                end;
            else
                % one-way ANOVA (paired) this is equivalent to unpaired t-test
                % -------------
                tail = 'one';
                [ori_vals df] = anova1_cell_select( data, g.paired );
                if strcmpi(g.method, 'param')
                    pvals = 1-fcdf(ori_vals, df(1), df(2));
                else
                    if strcmpi(g.arraycomp, 'on')
                        try
                            myfprintf(verb,'...');                        
                            res = surrogdistrib( data, 'method', g.method, 'pairing', g.paired, 'naccu', g.naccu);
                            surrogval = anova1_cell_select( res, g.paired );
                        catch,
                            myfprintf(verb,'\nSuperfast array computation failed because of memory limitation, reverting to standard computing');
                            g.arraycomp = 'off';
                        end;
                    end;
                    if strcmpi(g.arraycomp, 'off')
                        [res precomp] = surrogdistrib( data, 'method', g.method, 'pairing', g.paired);
                        for index = 1:g.naccu
                            if mod(index, 10) == 0, myfprintf(verb,'%d ', index); end;
                            if mod(index, 100) == 0, myfprintf(verb,'\n'); end;

                            res = surrogdistrib( {}, 'precomp', precomp);
                            if myndims(data{1}) == 1
                            	 surrogval(index)     = anova1_cell_select( res, g.paired );
                            else surrogval(:,index)   = anova1_cell_select( res, g.paired );
                            end;
                        end;
                    end;
                end;
            end;
        else
            % two-way ANOVA (paired or unpaired)
            % ----------------------------------
            tail = 'one';
            [ ori_vals{1} ori_vals{2} ori_vals{3} df{1} df{2} df{3} ] = anova2_cell_select( data, g.paired );
            if strcmpi(g.method, 'param')
                pvals{1} = 1-fcdf(ori_vals{1}, df{1}(1), df{1}(2));
                pvals{2} = 1-fcdf(ori_vals{2}, df{2}(1), df{2}(2));
                pvals{3} = 1-fcdf(ori_vals{3}, df{3}(1), df{3}(2));
            else
                surrogval = { surrogval surrogval surrogval };
                dataori   = data;
                if strcmpi(g.arraycomp, 'on')
                    try
                        myfprintf(verb,'...');
                        res = surrogdistrib( data, 'method', g.method, 'pairing', g.paired, 'naccu', g.naccu);
                        [ surrogval{1} surrogval{2} surrogval{3} ] = anova2_cell_select( res, g.paired );
                    catch,
                        myfprintf(verb,'\nSuperfast array computation failed because of memory limitation, reverting to standard computing');
                        g.arraycomp = 'off';
                    end;
                end;
                if strcmpi(g.arraycomp, 'off')
                    [res precomp] = surrogdistrib( data, 'method', g.method, 'pairing', g.paired);
                    for index = 1:g.naccu
                        if mod(index, 10) == 0, myfprintf(verb,'%d ', index); end;
                        if mod(index, 100) == 0, myfprintf(verb,'\n'); end;

                        res = surrogdistrib( {}, 'precomp', precomp);
                        if myndims(data{1}) == 1
                         	 [ surrogval{1}(index)     surrogval{2}(index)     surrogval{3}(index)     ] = anova2_cell_select( res, g.paired );
                        else [ surrogval{1}(:,index)   surrogval{2}(:,index)   surrogval{3}(:,index)   ] = anova2_cell_select( res, g.paired );
                        end;
                    end;
                end;
            end;
        end;
        myfprintf(verb,'\n');
    else
        surrogval = g.surrog;
        ori_vals  = g.stats;
        df        = [];
    end;
    
    % compute p-values
    % ----------------
    if ~strcmpi(g.method, 'param')
        if iscell( surrogval )
            pvals{1} = stat_surrogate_pvals(surrogval{1}, ori_vals{1}, tail);
            pvals{2} = stat_surrogate_pvals(surrogval{2}, ori_vals{2}, tail);
            pvals{3} = stat_surrogate_pvals(surrogval{3}, ori_vals{3}, tail);
        else
            pvals = stat_surrogate_pvals(surrogval, ori_vals, tail);
        end;
        try, warning('on', 'MATLAB:divideByZero'); catch, end;
    end;

    [ ori_vals, pvals ] = reshape_results( nd, ori_vals, pvals);
    [ surrogval ]       = reshape_results( [nd g.naccu], surrogval);
    
    % confidence intervals
    % --------------------
    if ~isnan(g.alpha)
        outputstruct.ci = stat_surrogate_ci(surrogval, g.alpha, tail);
        if strcmpi(g.structoutput, 'off')
            disp('Warning: returning confidence interval requires an output structure');
        end;
        if iscell(pvals)
            for ind = 1:length(pvals)
                outputstruct.mask{ind} = pvals{ind} < g.alpha;
            end;
        else
            outputstruct.mask = pvals < g.alpha;
        end;
    end;
    
    % create a structure for outputing values
    % ---------------------------------------
    if strcmpi(g.structoutput, 'on')
        outputstruct.method = g.method;
        outputstruct.pval   = pvals;
        outputstruct.df     = df;
        outputstruct.surrog = surrogval;
        if length(data(:)) == 2
             outputstruct.t = ori_vals;
        else outputstruct.f = ori_vals;
        end;
        outputstruct.stat   = ori_vals;
        ori_vals = outputstruct;
    end;
       
% compute ANOVA 2-way
% -------------------
function [f1 f2 f3 df1 df2 df3] = anova2_cell_select( res, paired);
    if strcmpi(paired,'on')
        [f1 f2 f3 df1 df2 df3] = anova2rm_cell( res );
    else
        [f1 f2 f3 df1 df2 df3] = anova2_cell( res );
    end;
    
% compute ANOVA 1-way
% -------------------
function [f df] = anova1_cell_select( res, paired);
    if strcmpi(paired,'on')
        [f df] = anova1rm_cell( res );
    else
        [f df] = anova1_cell( res );
    end;

% compute t-test
% -------------------
function [t df] = ttest_cell_select( res, paired, homogenous);
    if strcmpi(paired,'on')
        [t df] = ttest_cell( res{1}, res{2});
    else
        [t df] = ttest2_cell( res{1}, res{2}, homogenous);
    end;

% function to compute the number of dimensions
% --------------------------------------------
function val = myndims(a)
    if ndims(a) > 2
        val = ndims(a);
    else
        if size(a,1) == 1,
            val = 2;
        elseif size(a,2) == 1,
            val = 1;
        else
            val = 2;
        end;
    end; 

% function for verbose messages
% -----------------------------
function myfprintf(verb, varargin)
    if verb
        fprintf(varargin{:});
    end;

% function to replace tcdf
% ------------------------
function p = mytcdf(x,v)

if length(v) == 1,
    v = repmat(v, size(x));
end;

x2 = x.^2;
inds1 = (v < x2);
inds2 = (v >= x2);
if any(inds1(:)), p(inds1) = betainc(v(inds1) ./ (v(inds1) + x2(inds1)), v(inds1)/2, 0.5, 'lower') / 2; end;
if any(inds2(:)), p(inds2) = betainc(x2(inds2) ./ (v(inds2) + x2(inds2)), 0.5, v(inds2)/2, 'upper') / 2; end;
inds = (x > 0); 
if any(inds)
    p(inds) = 1 - p(inds);
end;

inds = (v > 1e7);
if any(inds(:)), p(inds) = normcum(x(inds)); end;

p(x == 0) = 0.5;
if isempty(p)
    p = ones(size(x));
else
    p = reshape(p, size(x));
end;
function [p] = normcum(z)
p = 0.5 * erfc(-z ./ sqrt(2));

% reshape results
% ---------------
function varargout = reshape_results(nd, varargin)
    if length(varargin) > 1
        for index = 1:length(varargin)
            varargout{index} = reshape_results(nd, varargin{index});
        end;
    elseif iscell(varargin{1})
        for index = 1:length(varargin{1})
            varargout{1}{index} = reshape_results(nd, varargin{1}{index});
        end;
    else
        if ~isempty(varargin{1})
            if length(nd) == 1, nd = [ nd 1 ]; end;
            varargout{1} = reshape(varargin{1}, nd);
        else varargout{1} = [];
        end;
    end;    
