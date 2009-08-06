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
%   'mode'     = ['perm'|'bootstrap'|'param'] mode for computing the p-values:
%                 'param' = parametric testing (standard ANOVA or t-test); 
%                 'perm' = non-parametric testing using surrogate data
%                 'bootstrap' = non-parametric bootstrap 
%                  made by permuting the input data {default: 'param'}
%   'naccu'    = [integer] Number of surrogate data copies to use in 'perm' 
%                 or 'bootstrap' mode estimation (see above) {default: 200}.
%   'verbose'  = ['on'|'off'] print info on the command line {default: 'on'}.
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
% Revision 1.20  2009/05/31 02:22:10  arno
% Adding FDR and bootstrap to all STUDY functions
%
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
                                 'mode'    'string'    { 'param' 'perm' 'bootstrap' }  'param';
                                 'paired'  'string'    { 'on' 'off' }      'on' 
                                 'verbose' 'string'    { 'on' 'off' }      'on' }, 'statcond');
    if isstr(g), error(g); end;
    
    if strcmpi(g.verbose, 'on'), verb = 1; else verb = 0; end;
    if strcmp(g.mode, 'param' ) & exist('fcdf') ~= 2
      myfprintf(verb,['statcond(): parametric testing requires fcdf() \n' ...
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
        
    % test if data can be paired
    % --------------------------
    if length(unique(cellfun('size', data, ndims(data{1}) ))) > 1
        g.paired = 'off'; 
    end;
    myfprintf(verb,'%d x %d, ', size(data,1), size(data,2));
    if strcmpi(g.paired, 'on')
         myfprintf(verb,'paired data, ');
         pairflag = 1;
    else myfprintf(verb,'unpaired data, ');
         pairflag = 0;
    end;
    if size(data,1) == 1 & size(data,2) == 2
         myfprintf(verb,'computing T values\n');
    else myfprintf(verb,'computing F values\n');
    end;
    
    % bootstrap flag
    % --------------
    if strcmpi(g.mode, 'bootstrap'), bootflag = 1;
    else                             bootflag = 0;
    end;
    
    % concatenate all data arrays
    % ---------------------------
    [ datavals datalen datadims ] = concatdata( data );
    
    % output text
    % -----------
    if ~strcmpi(g.mode, 'param')
        if bootflag, myfprintf(verb,'Bootstraps (of %d):', g.naccu);
        else         myfprintf(verb,'Permutations (of %d):', g.naccu);
        end;
    end;

    if size(data,1) == 1, % only one row
        
        if size(data,2) == 2 & strcmpi(g.paired, 'on')
            
            % paired t-test (very fast)
            % -------------
            tail = 'both';
            cond1 = data{1,1};
            cond2 = data{1,2};
            [ori_vals df] = paired_ttest(cond1, cond2);
            if strcmpi(g.mode, 'param')
                pvals = tcdf(ori_vals, df)*2;
                tmppvals = reshape(pvals, prod(size(pvals)),1);
                inds     = find(tmppvals > 1); 
                tmppvals(inds) = 2-tmppvals(inds);
                pvals    = reshape(tmppvals, size(pvals));
                return;
            else
                cond1ori = cond1;
                cond2ori = cond2;
                for index = 1:g.naccu
                    
                    res = surrogate( datavals, datalen, datadims, bootflag, pairflag);
                    cond1 = res{1,1};
                    cond2 = res{1,2};
                    if mod(index, 10) == 0, myfprintf(verb,'%d ', index); end;
                    if mod(index, 100) == 0, myfprintf(verb,'\n'); end;
                    switch myndims(cond1)
                     case 1   , surrogval(index)     = paired_ttest( cond1, cond2);
                     case 2   , surrogval(:,index)   = paired_ttest( cond1, cond2);
                     otherwise, surrogval(:,:,index) = paired_ttest( cond1, cond2);
                    end;
                    
                end;
            end;
        elseif size(data,2) == 2 & strcmpi(g.paired, 'off')
            
            % unpaired t-test (very fast)
            % -------------
            tail = 'both';
            [datac datalim ] = concatdata(data);
            
            cond1 = data{1,1};
            cond2 = data{1,2};
            [ori_vals df] = unpaired_ttest(cond1, cond2);
            if strcmpi(g.mode, 'param')
                pvals = tcdf(ori_vals, df)*2;
                tmppvals = reshape(pvals, prod(size(pvals)),1);
                inds     = find(tmppvals > 1); 
                tmppvals(inds) = 2-tmppvals(inds);
                pvals    = reshape(tmppvals, size(pvals));
                return;
            else
                for index = 1:g.naccu
                    
                    res = surrogate( datavals, datalen, datadims, bootflag, pairflag);
                    cond1 = res{1,1};
                    cond2 = res{1,2};
                    if mod(index, 10) == 0 , myfprintf(verb,'%d ', index); end;
                    if mod(index, 100) == 0, myfprintf(verb,'\n'); end;
                    switch myndims(cond1)
                     case 1   , surrogval(index)     = unpaired_ttest( cond1, cond2);
                     case 2   , surrogval(:,index)   = unpaired_ttest( cond1, cond2);
                     otherwise, surrogval(:,:,index) = unpaired_ttest( cond1, cond2);
                    end;
                    
                end;
            end;
            
        else
            % one-way ANOVA (paired) this is equivalent to unpaired t-test
            % -------------
            tail = 'one';
            [ori_vals df] = anova1_cell( data );
            if strcmpi(g.mode, 'param')
                pvals = 1-fcdf(ori_vals, df(1), df(2)); return;
            else
                for index = 1:g.naccu
                    
                    if mod(index, 10) == 0, myfprintf(verb,'%d ', index); end;
                    if mod(index, 100) == 0, myfprintf(verb,'\n'); end;
                    
                    res = surrogate( datavals, datalen, datadims, bootflag, pairflag);                    
                    switch myndims(data{1})
                     case 1   , surrogval(index)     = anova1_cell( res );
                     case 2   , surrogval(:,index)   = anova1_cell( res );
                     otherwise, surrogval(:,:,index) = anova1_cell( res );
                    end;
                    
                end;
            end;
        end;
    else
        % two-way ANOVA (paired)
        % ----------------------
        tail = 'one';
        [ ori_vals{1} ori_vals{2} ori_vals{3} df{1} df{2} df{3} ] = anova2_cell( data );
        if strcmpi(g.mode, 'param')
            pvals{1} = 1-fcdf(ori_vals{1}, df{1}(1), df{1}(2));
            pvals{2} = 1-fcdf(ori_vals{2}, df{2}(1), df{2}(2));
            pvals{3} = 1-fcdf(ori_vals{3}, df{3}(1), df{3}(2));
            return;
        else
            surrogval = { surrogval surrogval surrogval };
            dataori   = data;
            for index = 1:g.naccu
                
                if mod(index, 10) == 0, myfprintf(verb,'%d ', index); end;
                if mod(index, 100) == 0, myfprintf(verb,'\n'); end;
                
                res = surrogate( datavals, datalen, datadims, bootflag, pairflag);                    
                switch myndims(data{1})
                 case 1   , [ surrogval{1}(index)     surrogval{2}(index)     surrogval{3}(index)     ] = anova2_cell( res );
                 case 2   , [ surrogval{1}(:,index)   surrogval{2}(:,index)   surrogval{3}(:,index)   ] = anova2_cell( res );
                 otherwise, [ surrogval{1}(:,:,index) surrogval{2}(:,:,index) surrogval{3}(:,:,index) ] = anova2_cell( res );
                end;
                
            end;
        end;
    end;
    myfprintf(verb,'\n');
    
    % compute p-values
    % ----------------
    if iscell( surrogval )
        pvals{1} = compute_pvals(surrogval{1}, ori_vals{1}, tail);
        pvals{2} = compute_pvals(surrogval{2}, ori_vals{2}, tail);
        pvals{3} = compute_pvals(surrogval{3}, ori_vals{3}, tail);
    else
        pvals = compute_pvals(surrogval, ori_vals, tail);
    end;
    
% compute p-values
% ----------------
function pvals = compute_pvals(surrog, oridat, tail)
    
    surrog        = sort(surrog, myndims(surrog)); % sort last dimension
    
    if myndims(surrog) == 1    
        surrog(end+1) = oridat;        
    elseif myndims(surrog) == 2
        surrog(:,end+1) = oridat;        
    elseif myndims(surrog) == 3
        surrog(:,:,end+1) = oridat;
    else
        surrog(:,:,:,end+1) = oridat;
    end;

    [tmp idx] = sort( surrog, myndims(surrog) );
    [tmp mx]  = max( idx,[], myndims(surrog));        
                
    len = size(surrog,  myndims(surrog) );
    pvals = 1-(mx-0.5)/len;
    if strcmpi(tail, 'both')
        pvals = min(pvals, 1-pvals);
        pvals = 2*pvals;
    end;    
        
function res = surrogate(dataconcat, lens, dims, bootstrapflag, pairedflag); % for increased speed only shuffle half the indices
       
    % recompute indices in set and target cell indices
    % ------------------------------------------------
    if bootstrapflag
        if pairedflag
             indswap  = mod( [1:lens(end)]+ ceil(rand(1,lens(end))*length(lens))*lens(2)-1, lens(end) )+1;
        else indswap  = ceil(rand(1,lens(end))*lens(end));
        end;
    else
        if pairedflag
            indswap  = [1:lens(end)];
            indswap  = reshape(indswap, [lens(2) length(lens)-1]);
            for i = 1:size(indswap,1) % shuffle each row
                [tmp idx] = sort(rand(1,size(indswap,2)));
                indswap(i,:) = indswap(i,idx);
            end;    
            indswap  = reshape(indswap, [1 lens(2)*(length(lens)-1)]);
        else
            oriindices = [1:lens(end)]; % just shuffle indices
            [tmp idx] = sort(rand(1,length(oriindices)));
            indswap   = oriindices(idx);
        end;
    end;
    
    res = {};
    for i = 1:length(lens)-1
        switch myndims(dataconcat)
            case 1, res{i} = dataconcat(indswap(lens(i)+1:lens(i+1)));
            case 2, res{i} = dataconcat(:,indswap(lens(i)+1:lens(i+1)));
            case 3, res{i} = dataconcat(:,:,indswap(lens(i)+1:lens(i+1)));
            case 1, res{i} = dataconcat(:,:,:,indswap(lens(i)+1:lens(i+1)));
        end;
    end;
    res = reshape(res, dims);
    
function [tval, df] = paired_ttest(a,b)
    
    tmpdiff = a-b;
    diff = mymean(tmpdiff,    myndims(a));
    sd   = mystd( tmpdiff,[], myndims(a));
    tval = diff./sd*sqrt(size(a, myndims(a)));
    df   = size(a, myndims(a))-1;
    
    % check values againg Matlab statistics toolbox
    %[h p ci stats] = ttest(a', b');
    % [ tval stats.tstat' ]     
    
function [tval, df] = unpaired_ttest(a,b) % assumes equal variances
    
    meana = mymean(a, myndims(a));
    meanb = mymean(b, myndims(b));
    sda   = mystd(a, [], myndims(a));
    sdb   = mystd(b, [], myndims(b));
    na    = size(a, myndims(a));
    nb    = size(b, myndims(b));
    sp    = sqrt(((na-1)*sda.^2+(nb-1)*sdb.^2)/(na+nb-2));
    tval  = (meana-meanb)./sp/sqrt(1/na+1/nb);
    df    = na+nb-2;
    
    % check values againg Matlab statistics toolbox
    % [h p ci stats] = ttest2(a', b');
    % [ tval stats.tstat' ]     
            
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
  
function res = mymean( data, varargin) % deal with complex numbers
    res = mean( data, varargin{:});
    if ~isreal(data)
        res = abs( res );
    end;

function res = mystd( data, varargin) % deal with complex numbers
    if ~isreal(data)
        res = std( abs(data), varargin{:});
    else
        res = std( data, varargin{:});
    end;

function myfprintf(verb, varargin)
    if verb
        fprintf(varargin{:});
    end;
      