% statcond()  - compare two or more data conditions statistically using 
%               standard parametric or nonparametric permutation-based ANOVA 
%               (1-way or 2-way) or t-test methods. Parametric testing uses 
%               fcdf() from the Matlab Statistical Toolbox. Use of up to 
%               4-D data matrices speeds processing.
% Usage:
%          >> [stats, df, pvals, surrog] = statcond( data, 'key','val'... );
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
%                   The test used depends on the size of the data array input: 
%                When the data cell array has 2 columns and the data are 
%                paired, a paired t-test is performed; when the data are 
%                unpaired, an unpaired t-test is performed. If 'data' 
%                has only one row (paired or unpaired) and more than 2 
%                columns, a one-way ANOVA is performed. If the data cell 
%                array contains several rows and columns (of paired or 
%                unpaired data), a two-way ANOVA is performed.
%
% Optional inputs:
%   'paired'   = ['on'|'off'] pair the data array {default: 'on' unless 
%                the last dimension of data array is of different lengths}.
%   'mode'     = ['perm'|'param'] mode for computing the p-values:
%                 'param' = parametric testing (standard ANOVA or t-test); 
%                 'perm' = non-parametric testing using surrogate data 
%                  made by permuting the input data {default: 'perm'}
%   'naccu'    = [integer] Number of surrogate data copies to use in 'perm' 
%                 mode estimation (see above) {default: 200}.
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
%                 with three elements: output(1) = column effects; 
%                 output(2) = row effects; output(3) = interactions
%                 between rows and columns.
% Examples:
%      >> a = { rand(1,10) rand(1,10)+0.5 }; % pseudo 'paired' data vectors
%         [t df pvals] = statcond(a);        % perform paired t-test
%           pvals =                  
%              5.2807e-04 % standard t-test probability value
%         % Note: for different rand() outputs, results will differ.
%
%         [t df pvals surog] = statcond(a, 'mode', 'perm', 'naccu', 2000); 
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
%           pvals{1} % a (3,4) matrix of p-values; effects across columns
%           pvals{2} % a (3,4) matrix of p-values; effects across rows 
%           pvals{3} % a (3,4) matrix of p-values; interaction effects
%                                      % across rows and columns
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 2005-
%         With rhanks to Robert Oostenveld for fruitful discussions 
%         and advice on this function.
%
% See also: anova1_cell(), anova2_cell(), fcdf()

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
% [F df pval] = statcond(a, 'mode', 'perm', 'naccu', 200); pval
% [h p t stat] = ttest( a{1}(1,:), a{2}(1,:)); p

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.19  2007/05/04 23:19:02  arno
% *** empty log message ***
%
% Revision 1.18  2007/05/04 23:18:18  arno
% same
%
% Revision 1.17  2007/05/04 23:17:02  arno
% compute correct pval for unpaired t-test
%
% Revision 1.16  2007/04/27 22:06:34  arno
% unpaired condition shuffling
%
% Revision 1.15  2007/04/26 22:53:17  arno
% stat for t-test
%
% Revision 1.14  2007/04/06 19:34:42  arno
% Matlab 6.5 compatibiltiy
%
% Revision 1.13  2006/12/30 01:16:55  scott
% help msg
% quit
%
% Revision 1.12  2006/11/22 19:03:28  arno
% number of dimensions fix
%
% Revision 1.11  2006/11/22 18:46:15  arno
% fix typo
%
% Revision 1.10  2006/11/07 23:11:16  arno
% implement umpaired t-test, add documentation
%
% Revision 1.9  2006/10/03 22:00:58  scott
% minor help msg edits -sm
%
% Revision 1.8  2006/05/11 15:48:07  arno
% now better detect unpaired data
%
% Revision 1.5  2005/12/11 03:06:40  scott
% help msg editing -sm
%
% Revision 1.4  2005/12/09 19:05:13  arno
% editing header
%
% Revision 1.3  2005/12/09 18:43:39  arno
% scott's edit
%
% Revision 1.2  2005/12/07 23:31:23  arno
% chaging help message - allowing 4-D permutation
%
% Revision 1.1  2005/12/07 17:26:09  arno
% Initial revision
%

function [ ori_vals, df, pvals, surrogval ] = statcond( data, varargin );
    
    if nargin < 1
        help statcond;
        return;
    end;
    
    g = finputcheck( varargin, { 'naccu'   'integer'   [1 Inf]             200;
                                 'mode'    'string'    { 'param' 'perm' }  'param';
                                 'paired'  'string'    { 'on' 'off' }      'on' }, 'statcond');
    if isstr(g), error(g); end;
    
    if strcmp(g.mode, 'param' ) & exist('fcdf') ~= 2
      fprintf(['statcond(): parametric testing requires fcdf() \n' ...
               '            from the Matlab StatsticaL Toolbox.\n' ...
               '            Running nonparametric permutation tests\n.']);
      g.mode = 'perm';
    end
    if size(data,2) == 1, data  = transpose(data); end; % cell array transpose
    
    tmpsize   = size(data{1});
    if ~strcmpi(g.mode, 'param')
         surrogval = zeros([ tmpsize(1:end-1) g.naccu ]);
    else surrogval = [];
    end;
    
    % cfg configuration for Fieldtrip
    % -------------------------------
    cfg.method      = 'montecarlo';
    cfg.correctm    = 'no';
    cfg.neighbours  = {};
    cfg.feedback    = 'no';
    cfg.ivar        = 1;
    cfg.numrandomization = g.naccu;
            
    % test if data can be paired
    % --------------------------
    if length(unique(cellfun('size', data, ndims(data{1}) ))) > 1
        g.paired = 'off'; 
    end;
    fprintf('%d x %d, ', size(data,1), size(data,2));
    if strcmpi(g.paired, 'on')
         fprintf('paired data, ');
    else fprintf('unpaired data, ');
    end;
    if size(data,1) == 1 & size(data,2) == 2
         fprintf('computing T values\n');
    else fprintf('computing F values\n');
    end;
    
    % output text
    % -----------
    if ~strcmpi(g.mode, 'param')
        fprintf('Accumulating (of %d):', g.naccu);
    end;

    if size(data,1) == 1, % only one row
        
        if size(data,2) == 2 & strcmpi(g.paired, 'on')
            
            % paired t-test (very fast)
            % -------------
            cfg.statistic   = 'depsamplesT';
            [newdata design1 design2 design3] = makefieldtripdata(data, 0);
            cfg.design      = [ design1; design3 ];
            cfg.uvar        = 2;
            stat            = freqstatistics(cfg, newdata{:});
            ori_vals        = stat.stat;
            df              = [];
            pvals           = stat.prob;
            return;
            
        elseif size(data,2) == 2 & strcmpi(g.paired, 'off')
            
            % paired t-test (very fast)
            % -------------
            cfg.statistic   = 'indepsamplesT';
            [newdata design1] = makefieldtripdata(data, 0);
            cfg.design      = design1;            
            stat            = freqstatistics(cfg, newdata{:});
            ori_vals        = stat.stat;
            df              = [];
            pvals           = stat.prob;
            return;
            
        elseif strcmpi(g.paired, 'on')
            
            % one-way ANOVA (paired) this is equivalent to unpaired t-test
            % -------------
            cfg.statistic   = 'depsamplesF';
            [newdata design1 design2 design3] = makefieldtripdata(data, 0);
            cfg.design      = [ design1; design3 ];
            cfg.uvar        = 2;
            stat            = freqstatistics(cfg, newdata{:});
            ori_vals        = stat.stat;
            df              = [];
            pvals           = stat.prob;
            return;
            
        else
            % one-way ANOVA (unpaired) 
            % -------------
            cfg.statistic   = 'indepsamplesF';
            [newdata design1] = makefieldtripdata(data, 0);
            cfg.design      = [ design1 ];
            stat            = freqstatistics(cfg, newdata{:});
            ori_vals        = stat.stat;
            df              = [];
            pvals           = stat.prob;
            return;
            
        end;
    else
        if strcmpi(g.paired, 'on')
            
            % two-way ANOVA (paired) 
            % -------------
            cfg.statistic   = 'anovan';
            [newdata design1 design2 design3] = makefieldtripdata(data, 0);
            cfg.design      = [ design1; design2; design3 ];
            cfg.effect      = 'X1*X2';
            cfg.ivar        = [1 2];
            cfg.uvar        = 3;
            stat            = freqstatistics(cfg, newdata{:});
            ori_vals        = stat.stat;
            df              = [];
            pvals           = stat.prob;
            return;
            
        else
            
            % two-way ANOVA (unpaired)
            % -------------
            cfg.statistic   = 'anovan';
            [newdata design1 design2] = makefieldtripdata(data, 0);
            cfg.design      = [ design1; design2 ];
            cfg.effect      = 'X1*X2';
            cfg.ivar        = [1 2];
            stat            = freqstatistics(cfg, newdata{:});
            ori_vals        = stat.stat;
            df              = [];
            pvals           = stat.prob;
            return;
            
        end;
    end;
                
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
  
function [newdata, design1, design2, design3] = makefieldtripdata(data, chanflag);
    
    newdata = {};
    for i = 1:length(data(:))
    
        switch myndims(data{1})
          case 1, 
            newdata{i}.powspctrm = data{i};
            
          case 2,
            if chanflag
                 newdata{i}.powspctrm = transpose(data{i});
            else newdata{i}.powspctrm = reshape(transpose(data{i}), size(data{i},2), 1, size(data{i},1));
            end;
            
          case 3,
            if chanflag
                 newdata{i}.powspctrm = permute(data{i}, [3 1 2]);
            else newdata{i}.powspctrm = permute(data{i}, [3 4 1 2]);
            end;
            
          case 4,
            if chanflag
                 newdata{i}.powspctrm = permute(data{i}, [4 1 2 3]);
            else newdata{i}.powspctrm = permute(data{i}, [3 5 1 2]);
            end;

        end;
        
        newdata{i}.dimord    = 'rpt_chan_freq_time';
        newdata{i}.label     = cell(1,size(newdata{i}.powspctrm,2));
        newdata{i}.label(:)  = { 'cz' };
        newdata{i}.freq      = [1:size(newdata{i}.powspctrm,3)];
        newdata{i}.time      = [1:size(newdata{i}.powspctrm,4)];

    end;
    
    design1 = [];
    design2 = [];
    design3 = [];
    data
    for i = 1:size(data,2)
        for j = 1:size(data,1)
            ij = j+(i-1)*size(data,1);
            design1 = [ design1 ones(1, size(newdata{i}.powspctrm,1))*i ];
            design2 = [ design2 ones(1, size(newdata{i}.powspctrm,1))*j ];
            design3 = [ design3 [1:size(newdata{i}.powspctrm,1)] ];
        end;
    end;
