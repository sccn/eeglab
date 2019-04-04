% ttest2_cell() - compute unpaired t-test. Allow fast computation of 
%                 multiple t-test using matrix manipulation.
%
% Usage:
%    >> [F df] = ttest2_cell( { a b } );
%    >> [F df] = ttest2_cell(a, b);
%    >> [F df] = ttest2_cell(a, b, 'inhomogenous');
%
% Inputs:
%   a,b       = data consisting of UNPAIRED arrays to be compared. The last 
%               dimension of the data array is used to compute the t-test.
%   'inhomogenous' = use computation for the degree of freedom using 
%                    inhomogenous variance. By default the computation of
%                    the degree of freedom is done with homogenous
%                    variances.
%
% Outputs:
%   T   - T-value
%   df  - degree of freedom (array)
%
% Example:
%   a = { rand(1,10) rand(1,10)+0.5 }
%   [T df] = ttest2_cell(a)
%   signif = 2*tcdf(-abs(T), df(1))
%
%   % for comparison, the same using the Matlab t-test function
%   [h p ci stats] = ttest2(a{1}', a{2}');
%   [ stats.tstat' p] 
%
%   % fast computation (fMRI scanner volume 100x100x100 and 10 control 
%   % subjects and 12 test subjects). The computation itself takes 0.5 
%   % seconds instead of  half an hour using the standard approach (1000000 
%   % loops and Matlab  t-test function)
%   a = rand(100,100,100,10); b = rand(100,100,100,10);
%   [F df] = ttest_cell({ a b });
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 2005
%         (thank you to G. Rousselet for providing the formula for 
%         inhomogenous variances).
%
% Reference:
%   Schaum's outlines in statistics (3rd edition). 1999. Mc Graw-Hill.
%   Howel, Statistical Methods for Psychology. 2009. Wadsworth Publishing.

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

function [tval, df] = ttest2_cell(a,b,c) % assumes equal variances
    
    if nargin < 1
        help ttest2_cell;
        return;
    end
    
    homogenous = 'homogenous';
    if nargin > 1 && ischar(b)
        homogenous = b;
    end
    if nargin > 2 && ischar(c)
        homogenous = c;
    end
    if iscell(a), 
        b = a{2}; 
        a = a{1}; 
    end
    if ~strcmpi(homogenous, 'inhomogenous') && ~strcmpi(homogenous, 'homogenous')
        error('Value for homogenous parameter can only be ''homogenous'' or ''inhomogenous''');
    end

    nd    = myndims(a);
    na    = size(a, nd);
    nb    = size(b, nd);
    meana = mymean(a, nd);
    meanb = mymean(b, nd);
    
    if strcmpi(homogenous, 'inhomogenous')
        % inhomogenous variance from Howel, 2009, "Statistical Methods for Psychology"
        % thank you to G. Rousselet for providing these formulas
        m  = meana - meanb;
        s1 = var(a,0,nd) ./ na;
        s2 = var(b,0,nd) ./ nb;
        se = sqrt(s1 + s2);
        sd = sqrt([s1.*na, s2.*nb]);
        tval = m ./ se;

        df = ((s1 + s2).^2) ./ ((s1.^2 ./ (na-1) + s2.^2 ./ (nb-1)));
    else
        sda   = mystd(a, [], nd);
        sdb   = mystd(b, [], nd);
        sp    = sqrt(((na-1)*sda.^2+(nb-1)*sdb.^2)/(na+nb-2));
        tval  = (meana-meanb)./sp/sqrt(1/na+1/nb);
        df    = na+nb-2;
    end
        
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
        end
    end; 
  
function res = mymean( data, varargin) % deal with complex numbers
    res = mean( data, varargin{:});
    if ~isreal(data)
        res = abs( res );
    end

function res = mystd( data, varargin) % deal with complex numbers
    if ~isreal(data)
        res = std( abs(data), varargin{:});
    else
        res = sqrt(sum( bsxfun(@minus, data, mean( data, varargin{2})).^2, varargin{2})/(size(data,varargin{2})-1)); % 8 percent speedup
        %res = std( data, varargin{:});
    end
    
    
