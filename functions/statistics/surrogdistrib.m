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
  
function [res precomp ] = surrogdistrib(data, varargin)

if nargin < 1
    help surrogdistrib,
    return;
end;

if ~strcmpi(varargin{1}, 'precomp')
    opt = finputcheck(varargin, { 'naccu'      'integer'   [1 Inf]             1;
                                  'method'     'string'    { 'perm','permutation','bootstrap' }  'perm';
                                  'pairing'    'string'    { 'on','off' }      'on';
                                  'precomp'    'cell'      {} {} }, 'surrogdistrib');
    if isstr(opt), error(opt); end;
    if strcmpi(opt.method, 'permutation'), opt.method = 'perm'; end;
    if strcmpi(opt.method, 'bootstrap'), bootflag = 1;
    else                                 bootflag = 0;
    end;
    if strcmpi(opt.pairing, 'on')
         pairflag = 1;
    else pairflag = 0;
    end;
else
    opt.precomp = varargin{2};
end;

% concatenate data
% ----------------
if isempty(opt.precomp)
    [ datavals datalen datadims ] = concatdata( data );
    precomp = { datavals datalen datadims bootflag pairflag opt.naccu};
else
    precomp   = opt.precomp;
    datavals  = precomp{1};
    datalen   = precomp{2};
    datadims  = precomp{3};
    bootflag  = precomp{4};
    pairflag  = precomp{5};
    opt.naccu = precomp{6};
end;

% compute surrogate distribution
% ------------------------------
if opt.naccu > 1
    res = supersurrogate( datavals, datalen, datadims, bootflag, pairflag, opt.naccu);
else
    res = surrogate( datavals, datalen, datadims, bootflag, pairflag);
end;

function res = supersurrogate(dat, lens, dims, bootstrapflag, pairedflag, naccu); % for increased speed only shuffle half the indices

    % recompute indices in set and target cell indices
    % ------------------------------------------------
    ncond = length(lens)-1;
    nsubj = lens(2);
    if bootstrapflag
        if pairedflag
             indswap  = mod( repmat([1:lens(end)],[naccu 1]) + ceil(rand(naccu,lens(end))*length(lens))*lens(2)-1, lens(end) )+1;
        else indswap  = ceil(rand(naccu,lens(end))*lens(end));
        end;
    else
        if pairedflag
            [tmp idx] = sort(rand(naccu,nsubj,ncond),3);
            indswap   = ((idx)-1)*nsubj + repmat( repmat([1:nsubj], [naccu 1 1]),[1 1 ncond]);
            indswap   = reshape(indswap, [naccu lens(end)]);
        else
            [tmp indswap] = sort(rand(naccu, lens(end)),2);
        end;
    end;

    for i = 1:length(lens)-1
        if myndims(dat) == 1
             res{i} = reshape(dat(indswap(:,lens(i)+1:lens(i+1))), naccu, lens(i+1)-lens(i));
        else res{i} = reshape(dat(:,indswap(:,lens(i)+1:lens(i+1))), size(dat,1), naccu, lens(i+1)-lens(i));
        end;
    end;
    res = reshape(res, dims);

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
        if myndims(dataconcat) == 1
             res{i} = dataconcat(indswap(lens(i)+1:lens(i+1)));
        else res{i} = dataconcat(:,indswap(lens(i)+1:lens(i+1)));
        end;
    end;
    res = reshape(res, dims);
    
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
    