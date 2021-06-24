% surrogdistrib - Build surrogate distribution
%
% surrog = surrogdistrib(data, varargin);
%
% Inputs:
%  data   - [cell] data arrays for which to compute a surrogate
%           distribution.
%
% Optional inputs:
%  'method'  - ['bootstrap'|'perm'] use either 'bootstrap' or 'permutation'
%              method. Bootstrap performs draws with replacement and 
%              permutation performs draws without replacement. Default 
%              is 'perm'. 
%  'pairing' - ['on'|'off'] pair the data arrays.
%  'naccu'   - [integer] number of surrogate. Default is 1.
%  'precomp' - cell array containing precomputed value for speeding up
%              mulitple calls
%
% Output:
%  surrog  - surrogate distribution
%  precomp - cell array containing precomputed value for speeding up
%            mulitple calls
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 2005-

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
  
function [res, precomp ] = surrogdistrib(data, varargin)

if nargin < 1
    help surrogdistrib,
    return;
end

if ~strcmpi(varargin{1}, 'precomp')
    opt = finputcheck(varargin, { 'naccu'      'integer'   [1 Inf]             1;
                                  'method'     'string'    { 'perm','permutation','bootstrap' }  'perm';
                                  'pairing'    'string'    { 'on' 'off' }      'on';
                                  'precomp'    'cell'      {} {} }, 'surrogdistrib');
    if ischar(opt), error(opt); end
    if strcmpi(opt.method, 'permutation'), opt.method = 'perm'; end
    if strcmpi(opt.method, 'bootstrap'), bootflag = 1;
    else                                 bootflag = 0;
    end
    if strcmpi(opt.pairing, 'on')
         pairflag = 1;
    else pairflag = 0;
    end
else
    opt.precomp = varargin{2};
end

% concatenate data
% ----------------
if isempty(opt.precomp)
    [ datavals, datalen, datadims ] = concatdata( data );
    precomp = { datavals datalen datadims bootflag pairflag opt.naccu};
else
    precomp   = opt.precomp;
    datavals  = precomp{1};
    datalen   = precomp{2};
    datadims  = precomp{3};
    bootflag  = precomp{4};
    pairflag  = precomp{5};
    opt.naccu = precomp{6};
end

% compute surrogate distribution
% ------------------------------
if opt.naccu > 1
    res = supersurrogate( datavals, datalen, datadims, bootflag, pairflag, opt.naccu);
else
    res = surrogate( datavals, datalen, datadims, bootflag, pairflag);
end

function res = supersurrogate(dat, lens, dims, bootstrapflag, pairedflag, naccu) % for increased speed only shuffle half the indices

    % recompute indices in set and target cell indices
    % ------------------------------------------------
    ncond = length(lens)-1;
    nsubj = lens(2);
    if bootstrapflag
        if pairedflag
             indswap  = mod( repmat([1:lens(end)],[naccu 1]) + ceil(rand(naccu,lens(end))*length(lens))*lens(2)-1, lens(end) )+1;
        else indswap  = ceil(rand(naccu,lens(end))*lens(end));
        end
    else
        if pairedflag
            [tmp, idx] = sort(rand(naccu,nsubj,ncond),3);
            indswap   = ((idx)-1)*nsubj + repmat( repmat([1:nsubj], [naccu 1 1]),[1 1 ncond]);
            indswap   = reshape(indswap, [naccu lens(end)]);
        else
            [tmp, indswap] = sort(rand(naccu, lens(end)),2);
        end
    end

    for i = 1:length(lens)-1
        if myndims(dat) == 1
             res{i} = reshape(dat(indswap(:,lens(i)+1:lens(i+1))), naccu, lens(i+1)-lens(i));
        else res{i} = reshape(dat(:,indswap(:,lens(i)+1:lens(i+1))), size(dat,1), naccu, lens(i+1)-lens(i));
        end
    end
    res = reshape(res, dims);

function res = surrogate(dataconcat, lens, dims, bootstrapflag, pairedflag) % for increased speed only shuffle half the indices
       
    % recompute indices in set and target cell indices
    % ------------------------------------------------
    if bootstrapflag
        if pairedflag
             indswap  = mod( [1:lens(end)]+ ceil(rand(1,lens(end))*length(lens))*lens(2)-1, lens(end) )+1;
        else indswap  = ceil(rand(1,lens(end))*lens(end));
        end
    else
        if pairedflag
            indswap  = [1:lens(end)];
            indswap  = reshape(indswap, [lens(2) length(lens)-1]);
            for i = 1:size(indswap,1) % shuffle each row
                [tmp, idx] = sort(rand(1,size(indswap,2)));
                indswap(i,:) = indswap(i,idx);
            end;    
            indswap  = reshape(indswap, [1 lens(2)*(length(lens)-1)]);
        else
            oriindices = [1:lens(end)]; % just shuffle indices
            [tmp, idx] = sort(rand(1,length(oriindices)));
            indswap   = oriindices(idx);
        end
    end
    
    res = {};
    for i = 1:length(lens)-1
        if myndims(dataconcat) == 1
             res{i} = dataconcat(indswap(lens(i)+1:lens(i+1)));
        else res{i} = dataconcat(:,indswap(lens(i)+1:lens(i+1)));
        end
    end
    res = reshape(res, dims);
    
function val = myndims(a)
    if ndims(a) > 2
        val = ndims(a);
    else
        if size(a,1) == 1
            val = 2;
        elseif size(a,2) == 1
            val = 1;
        else
            val = 2;
        end
    end
    
