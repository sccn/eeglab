%  statcond() - compare conditions using surrogate t-test or ANOVA.
%
% Usage:
%    >> statcond( data, 'key', 'val', ... );
%
% Inputs:
%   data       = data consisting of arrays to be compared. The last dimension
%                of the data array is shufled across conditions.
%
% Optional inputs:
%   'paired'   = ['on'|'off'] pair the data array. Default is 'on'. 'off'
%                option not implemented yet.
%   'mode'     = ['perm'|'param'] mode for the computation of the p-value. 
%                'param is parametric (standard ANOVA) and 'perm' is 
%                using permutations. Default is 'perm'.
%   'naccu'    = [integer] for the 'perm' mode above, use that many 
%                repetition. Default is 200.
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 2005
%
% See also: statcond()

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

function [ ori_vals, df, pvals, surrogval ] = statcond( data, varargin );
    
    if nargin < 1
        help statcond;
        return;
    end;
    
    g = finputcheck( varargin, { 'naccu'   'integer'   [1 Inf]             200;
                                 'mode'    'string'    { 'param' 'perm' }  'param';
                                 'paired'  'string'    { 'on' 'off' }      'on' }, 'statcond');
    if isstr(g), error(g); end;
    if strcmpi(g.paired, 'off'), error('Paired ANOVA not implemented yet'); end;
    
    if size(data,2) == 1, data  = transpose(data); end;
    
    tmpsize   = size(data{1});
    surrogval = zeros([ tmpsize(1:end-1) g.naccu ]);
        
    % test if data can be paired
    % --------------------------
    if length(unique(cellfun('length', data ))) > 1
        g.paired = 'off'; 
    end;
    if strcmpi(g.paired, 'on')
         disp('Paired data');
    else disp('Unpaired data');
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
        else
            tmp            = a(:,:,indswap);
            a(:,:,indswap) = b(:,:,indswap);
            b(:,:,indswap) = tmp;
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
            else
                tmp                           = a{indarg1(i)}(:,:,indswap(i));
                a{indarg1(i)}(:,:,indswap(i)) = a{indarg2(i)}(:,:,indswap(i));
                a{indarg2(i)}(:,:,indswap(i)) = tmp;
            end;
        end;
        a = reshape(a, dims);
        
    end;
    

 function a = shuffle_unpaired(a); % for increased speed only shuffle half the indices
       
    % unpaired
    % --------
    dims      = size(a);
    a         = a(:)';                
    alllen       = cellfun('length', a );
    
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
        else
            tmp                             = a{indtarg1(i)}(:,:,indswap1(i));
            a{indtarg1(i)}(:,:,indswap1(i)) = a{indtarg2(i)}(:,:,indswap2(i));
            a{indtarg2(i)}(:,:,indswap2(i)) = tmp;
        end;
    end;
    a = reshape(a, dims);
        
function [tval, df] = paired_ttest(a,b)
    
    tmpdiff = a-b;
    diff = mean(tmpdiff,    myndims(a));
    sd   = std( tmpdiff,[], myndims(a));
    tval = diff./sd*sqrt(size(a, myndims(a)));
    df   = size(a, myndims(a))-1;
            
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
  
