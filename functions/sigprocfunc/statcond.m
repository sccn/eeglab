%  statcond()  - compare two or more data conditions statistically using standard 
%                parametric or nonparametric permutation-based ANOVA 
%                (1-way or 2-way) or t-test methods. Parametric testing uses 
%                fcdf() from the Matlab Statistical Toolbox. Use of up to 
%                4-D data matrices speeds processing.
% Usage:
%             >> [stats, df, pvals, surrog] = statcond( data, 'key', 'val', ... );
% Inputs:
%   data       = one-dim or two-dim cell array of data matrices. 
%                For nonparametric, permutation-based testing, the last 
%                dimension of the data arrays (which may be of up to 4
%                dimensions) is permuted across conditions, either in
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
%                the last dimension. A T value is returned only if the data 
%                consist of two conditions.
%   df         = degrees of freedom, a (2,1) vector if F-values are returned
%   pvals      = array of p-values. Same size as input data without the last
%                data dimension. All returned p-values are two tailed.
%   surrog     = surrogate data array (same size as input data with the last 
%                dim. filled with a number ('naccu') of surrogate data sets.
%
% Important note: When a two-way ANOVA is performed, outputs are cell arrays
%                 with three elements: output(1) = column effects; 
%                 output(2) = row effects; output(3) = interactions
%                 between rows and columns.
% Examples:
%      >> a = { rand(1,10) rand(1,10)+0.5 }; % pseudo 'paired' data vectors
%         [t df pvals] = statcond(a);        % perform paired t-test
%         pvals =                  
%            5.2807e-04          % standard t-test probability value
%            % Note: for different rand() outputs, results will differ.
%         [t df pvals surog] = statcond(a, 'mode', 'perm', 'naccu', 2000); 
%         pvals =
%            0.0065 % nonparametric t-test using 2000 permuted data sets
%
%         a = { rand(2,11) rand(2,10) rand(2,12)+0.5 }; % pseudo 'unpaired' 
%         [F df pvals] = statcond(a); % perform an unpaired ANOVA 
%         pvals =
%            0.00025 % p-values for difference between columns 
%            0.00002 % for each data row
%
%         a = { rand(3,4,10) rand(3,4,10) rand(3,4,10); ...
%               rand(3,4,10) rand(3,4,10) rand(3,4,10)+0.5 }; 
%         % pseudo (2,3)-condition data array, each entry containing 
%         %                                    ten (3,4) data matrices
%         [F df pvals] = statcond(a);  % perform a paired 2-way ANOVA 
%         pvals{1} % a (3,4) matrix of p-values; effects across columns
%         pvals{2} % a (3,4) matrix of p-values; effects across rows 
%         pvals{3} % a (3,4) matrix of p-values; interaction effects
%                                             % across rows and columns
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 2005-
%         Thanks to Robert Oostenveld for fruitful discussions and advices
%         about this function.
%
% See also: anova1_cell(), anova2_cell(), fcdf()

% testing paired t-test
% a = { rand(2,10) rand(2,10) };
% [t df pval] = statcond(a); pval
% [h p t stat] = ttest( a{1}(1,:), a{2}(1,:)); p
% [h p t stat] = ttest( a{1}(2,:), a{2}(2,:)); p
%
% compare significance level
% --------------------------
% a = { rand(1,10) rand(1,10) }; [F df pval] = statcond(a, 'mode', 'perm', 'naccu', 200); pval
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
            tail = 'both';
            cond1 = data{1,1};
            cond2 = data{1,2};
            [ori_vals df] = paired_ttest(cond1, cond2);
            if strcmpi(g.mode, 'param')
                pvals = tcdf(ori_vals, df)*2; return;
            else
                for index = 1:g.naccu
                    
                    [cond1 cond2]    = shuffle_paired(cond1, cond2);
                    if mod(index, 10) == 0, fprintf('%d ', index); end;
                    if mod(index, 100) == 0, fprintf('\n'); end;
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
            cond1 = data{1,1};
            cond2 = data{1,2};
            [ori_vals df] = unpaired_ttest(cond1, cond2);
            if strcmpi(g.mode, 'param')
                pvals = tcdf(ori_vals, df)*2; return;
            else
                for index = 1:g.naccu
                    
                    [cond1 cond2]    = shuffle_unpaired(cond1, cond2);
                    if mod(index, 10) == 0 , fprintf('%d ', index); end;
                    if mod(index, 100) == 0, fprintf('\n'); end;
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
                    
                    if mod(index, 10) == 0, fprintf('%d ', index); end;
                    if mod(index, 100) == 0, fprintf('\n'); end;
                    
                    if strcmpi(g.paired, 'on')
                         data = shuffle_paired(   data );
                    else data = shuffle_unpaired( data );
                    end;
                    
                    switch myndims(data{1})
                     case 1   , surrogval(index)     = anova1_cell( data );
                     case 2   , surrogval(:,index)   = anova1_cell( data );
                     otherwise, surrogval(:,:,index) = anova1_cell( data );
                    end;
                    
                end;
            end;
        end;
    else
        % two-way ANOVA (paired)
        % -------------
        tail = 'one';
        [ ori_vals{1} ori_vals{2} ori_vals{3} df{1} df{2} df{3} ] = anova2_cell( data );
        if strcmpi(g.mode, 'param')
            pvals{1} = 1-fcdf(ori_vals{1}, df{1}(1), df{1}(2));
            pvals{2} = 1-fcdf(ori_vals{2}, df{2}(1), df{2}(2));
            pvals{3} = 1-fcdf(ori_vals{3}, df{3}(1), df{3}(2));
            return;
        else
            surrogval = { surrogval surrogval surrogval };
            for index = 1:g.naccu
                
                if mod(index, 10) == 0, fprintf('%d ', index); end;
                if mod(index, 100) == 0, fprintf('\n'); end;
                
                if strcmpi(g.paired, 'on')
                     data = shuffle_paired(   data );
                else data = shuffle_unpaired( data );
                end;
                
                switch myndims(data{1})
                 case 1   , [ surrogval{1}(index)     surrogval{2}(index)     surrogval{3}(index)     ] = anova2_cell( data );
                 case 2   , [ surrogval{1}(:,index)   surrogval{2}(:,index)   surrogval{3}(:,index)   ] = anova2_cell( data );
                 otherwise, [ surrogval{1}(:,:,index) surrogval{2}(:,:,index) surrogval{3}(:,:,index) ] = anova2_cell( data );
                end;
                
            end;
        end;
    end;
    fprintf('\n');
    
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
    
% shuffle last dimension of arrays
% --------------------------------
function [a b] = shuffle_paired(a, b); % for increased speed only shuffle half the indices
    

    if nargin > 1 % very fast 2 conditions -> 2 inputs and 2 outputs
        
        indswap = find(round(rand( 1,size(a,myndims(a)) )));
        if myndims(a) == 1
            tmp        = a(indswap);
            a(indswap) = b(indswap);
            b(indswap) = tmp;
        elseif myndims(a) == 2    
            tmp          = a(:,indswap);
            a(:,indswap) = b(:,indswap);
            b(:,indswap) = tmp;
        elseif myndims(a) == 3
            tmp            = a(:,:,indswap);
            a(:,:,indswap) = b(:,:,indswap);
            b(:,:,indswap) = tmp;
        else
            tmp              = a(:,:,:,indswap);
            a(:,:,:,indswap) = b(:,:,:,indswap);
            b(:,:,:,indswap) = tmp;
        end;
        
    else % more than 2 conditions -> one cell array input and one cell array output
        
        dims      = size(a);
        a         = a(:)';
        in1       = a{1};
        indswap   = rand( length(a)-1,size(in1,myndims(in1)) ); % array of 0 and 1
        [tmp1 indswap] = find( indswap );                       % only shuffle half of them
        indarg1        = ceil(rand(1,length(indswap))*length(a));
        indarg2        = ceil(rand(1,length(indswap))*length(a));
    
        for i = 1:length(indswap)
            if myndims(in1) == 1
                tmp                       = a{indarg1(i)}(indswap(i)); % do not shuffle index
                a{indarg1(i)}(indswap(i)) = a{indarg2(i)}(indswap(i));
                a{indarg2(i)}(indswap(i)) = tmp;
            elseif myndims(in1) == 2    
                tmp                         = a{indarg1(i)}(:,indswap(i));
                a{indarg1(i)}(:,indswap(i)) = a{indarg2(i)}(:,indswap(i));
                a{indarg2(i)}(:,indswap(i)) = tmp;
            elseif myndims(in1) == 3
                tmp                           = a{indarg1(i)}(:,:,indswap(i));
                a{indarg1(i)}(:,:,indswap(i)) = a{indarg2(i)}(:,:,indswap(i));
                a{indarg2(i)}(:,:,indswap(i)) = tmp;
            else
                tmp                             = a{indarg1(i)}(:,:,:,indswap(i));
                a{indarg1(i)}(:,:,:,indswap(i)) = a{indarg2(i)}(:,:,:,indswap(i));
                a{indarg2(i)}(:,:,:,indswap(i)) = tmp;
            end;
        end;
        a = reshape(a, dims);
        
    end;
    

 function a = shuffle_unpaired(a); % for increased speed only shuffle half the indices
       
    % unpaired
    % --------
    dims      = size(a);
    a         = a(:)';                
    alllen    = cellfun('size', a, myndims(a{1}) ); % by chance, pick up the last dimension
    
    indswap1 = ceil(rand( 1, ceil(sum(alllen)/2) )*sum(alllen)); % origin index (cumulated indices)
    indswap2 = ceil(rand( 1, ceil(sum(alllen)/2) )*sum(alllen)); % origin target
    indtarg1 = ones(1, sum(alllen) );                    % origin set
    indtarg2 = ones(1, sum(alllen) );                    % target set
    
    % recompute indices in set and target cell indices
    % ------------------------------------------------
    for index = 1:(length(a)-1)
        tmpind1 = find(indswap1 > alllen(index));
        tmpind2 = find(indswap2 > alllen(index));
        indswap1(tmpind1) = indswap1(tmpind1)-alllen(index);
        indswap2(tmpind2) = indswap2(tmpind2)-alllen(index);
        indtarg1(tmpind1) = indtarg1(tmpind1)+1;
        indtarg2(tmpind2) = indtarg2(tmpind2)+1;
    end;
    
    % perform swaping
    % ---------------
    for i = 1:length(indswap1)
        if myndims(a{1}) == 1
            tmp                         = a{indtarg1(i)}(indswap1(i));
            a{indtarg1(i)}(indswap1(i)) = a{indtarg2(i)}(indswap2(i));
            a{indtarg2(i)}(indswap2(i)) = tmp;
        elseif myndims(a{1}) == 2    
            tmp                           = a{indtarg1(i)}(:,indswap1(i));
            a{indtarg1(i)}(:,indswap1(i)) = a{indtarg2(i)}(:,indswap2(i));
            a{indtarg2(i)}(:,indswap2(i)) = tmp;
        elseif myndims(a{1}) == 3    
            tmp                             = a{indtarg1(i)}(:,:,indswap1(i));
            a{indtarg1(i)}(:,:,indswap1(i)) = a{indtarg2(i)}(:,:,indswap2(i));
            a{indtarg2(i)}(:,:,indswap2(i)) = tmp;
        else
            tmp                               = a{indtarg1(i)}(:,:,:,indswap1(i));
            a{indtarg1(i)}(:,:,:,indswap1(i)) = a{indtarg2(i)}(:,:,:,indswap2(i));
            a{indtarg2(i)}(:,:,:,indswap2(i)) = tmp;
        end;
    end;
    a = reshape(a, dims);
        
function [tval, df] = paired_ttest(a,b)
    
    tmpdiff = a-b;
    diff = mymean(tmpdiff,    myndims(a));
    sd   = mystd( tmpdiff,[], myndims(a));
    tval = diff./sd*sqrt(size(a, myndims(a)));
    df   = size(a, myndims(a))-1;
            
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
    res = std( abs(data), varargin{:});
