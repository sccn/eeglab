% anova1_cell() - compute F-values in cell array using ANOVA.
%
% Usage:
%    >> [F df] = anova1_cell( data );
%
% Inputs:
%   data       = data consisting of PAIRED arrays to be compared. The last 
%                dimension of the data array is used to compute ANOVA.
% Outputs:
%   F   - F-value
%   df  - degree of freedom (array)
%
% Note: the advantage over the ANOVA1 function of Matlab statistical
%       toolbox is that this function works on arrays (see examples). Note
%       also that you still need the statistical toolbox to assess
%       significance using the fcdf() function. The other advantage is that
%       this function will work with complex numbers.
%
% Example:
%   a = { rand(1,10) rand(1,10) rand(1,10) }
%   [F df] = anova1_cell(a)
%   signif = 1-fcdf(F, df(1), df(2))
%
%   % for comparison 
%   anova1( [ a{1,1}' a{1,2}' a{1,3}' ]) % look in the graph for the F value
%
%   b = { [ a{1,1}; a{1,1} ] [ a{1,2}; a{1,2} ] [ a{1,3}; a{1,3} ] }
%   [F df] = anova1_cell(b)
%
%   c{1,1} = reshape(repmat(b{1,1}, [2 1]),2,2,10);
%   c{1,2} = reshape(repmat(b{1,2}, [2 1]),2,2,10);
%   c{1,3} = reshape(repmat(b{1,3}, [2 1]),2,2,10);
%   [F df] = anova1_cell(c)
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 2005
%
% Reference:
%   Schaum's outlines in statistics (3rd edition). 1999. Mc Graw-Hill.

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

function [F, df] = anova1_cell(data)
    
    % This function does not return
    % correct values (see bug 336)
    % It should be fixed with Schaum's outlines p363
    % but requires some work. It now calls
    % anova2_cell which returns correct values
    
    warning off;
    [ F tmp tmp2 df] =  anova2_cell(data);
    warning on;
    return;
    
    % compute all means and all std
    % -----------------------------
    nd = myndims( data{1} );
    if nd == 1
        
        for i = 1:length(data)
            n( i) = length(data{i});
            m( i) = mymean(  data{i});
            sd(i) = mystd(   data{i});
        end
        nt = sum(n);
        n   = n';
        m   = m';
        sd  = sd';
        
    elseif nd == 2  

        for i = 1:length(data)
            n( :,i) = ones(size(data{i},1) * size(data{i},2), 'single');
            m( :,i) = mymean(  data{i},2);
            sd(:,i) = mystd(   data{i},[],2);
        end
        nt = sum(n(1,:));
    
    elseif nd == 3        
        
        for i = 1:length(data)
            n( :,:,i) = ones(size(data{i},1),size(data{i},2) * size(data{i},3), 'single');
            m( :,:,i) = mymean(  data{i},3);
            sd(:,:,i) = mystd(   data{i},[],3);
        end
        nt = sum(n(1,1,:));
        
    elseif nd == 4

        for i = 1:length(data)
            n( :,:,:,i) = ones(size(data{i},1),size(data{i},2), size(data{i},3) * size(data{i},4), 'single');
            m( :,:,:,i) = mymean(  data{i},4);
            sd(:,:,:,i) = mystd(   data{i},[],4);
        end
        nt = sum(n(1,1,1,:));
        
    end
    
    mt = mean(m,nd);
    ng = length(data); % number of conditions
    
    VinterG  = ( sum( n.*(m.^2), nd ) - nt*mt.^2 )/(ng-1);
    VwithinG = sum( (n-1).*(sd.^2), nd )/(nt-ng);
    F  = VinterG./VwithinG;
    df = [ ng-1 ng*(size(data{1},nd)-1) ];

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
    end

function res = mymean( data, varargin) % deal with complex numbers
    res = mean( data, varargin{:});
    if ~isreal(data)
        res = abs( res );
    end

function res = mystd( data, varargin) % deal with complex numbers
    res = std( abs(data), varargin{:});
