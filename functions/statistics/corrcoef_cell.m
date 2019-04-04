% corrcoef_cell() - compute pairwise correlations using arrays and 
%                   cell array inputs.
%
% Usage:
%    >> c = corrcoef_cell( data );
%    >> c = corrcoef_cell( data );
%
% Inputs:
%   data       - [cell array] data consisting of PAIRED arrays to be compared. 
%                The last dimension of embeded data arrays is used to compute 
%                correlation (see examples).
% Outputs:
%   c   - Correlation values. Same size as data without the last dimension.
%
% Note: the main advantage over the corrcoef Matlab function is the
%       capacity to compute millions of pairwise correlations per second.
%
% Example:
%   a = { rand(1,10) rand(1,10) };
%   c1 = corrcoef_cell(a);
%   c2 = corrcoef(a{1}, a{2});
%   % in this case, c1 is equal to c2(2)
%
%   a = { rand(200,300,100) rand(200,300,100) };
%   c = corrcoef_cell(a);
%   % the call above would require 200 x 300 calls to the corrcoef function
%   % and be about 1000 times slower
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 2010

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

function c = corrcoef_cell(a,b);

if nargin < 1
    help corrcoef_cell;
    return;
end

if nargin < 2
    b = a{2};
    a = a{1};
end

nd = myndims(a);
if nd == 1
    aa = a-mean(a);
    bb = b-mean(b);
    cv  = aa'*bb/(10-1);
    cva = aa'*aa/(10-1);
    cvb = bb'*bb/(10-1);

    c = cv/sqrt(cva*cvb);
elseif nd == 2 % ND=2, 3, and 4 could be replaced with a single line
    aa = bsxfun(@minus, a, mean(a,2));
    bb = bsxfun(@minus, b, mean(b,2));
    %aa = a-repmat(mean(a,2),[1 size(a,2)]);
    %bb = b-repmat(mean(b,2),[1 size(a,2)]);
    cv  = sum(aa.*bb,2);
    cva = sum(aa.*aa,2);
    cvb = sum(bb.*bb,2);

    c = cv./sqrt(cva.*cvb);
elseif nd == 3
    aa = bsxfun(@minus, a, mean(a,3));
    bb = bsxfun(@minus, b, mean(b,3));
    %aa = a-repmat(mean(a,3),[1 1 size(a,3)]);
    %bb = b-repmat(mean(b,3),[1 1 size(a,3)]);
    cv  = sum(aa.*bb,3);
    cva = sum(aa.*aa,3);
    cvb = sum(bb.*bb,3);

    c = cv./sqrt(cva.*cvb);
elseif nd == 4
    aa = bsxfun(@minus, a, mean(a,4));
    bb = bsxfun(@minus, b, mean(b,4));
    %aa = a-repmat(mean(a,4),[1 1 1 size(a,4)]);
    %bb = b-repmat(mean(b,4),[1 1 1 size(a,4)]);
    cv  = sum(aa.*bb,4);
    cva = sum(aa.*aa,4);
    cvb = sum(bb.*bb,4);

    c = cv./sqrt(cva.*cvb);
end

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
