%  anova1rm_cell() - compute F-values in cell array using repeated measure
%                    ANOVA.
%
% Usage:
%    >> [FC dfc] = anova2rm_cell( data );
%
% Inputs:
%   data       = data consisting of PAIRED arrays to be compared. The last 
%                dimension of the data array is used to compute ANOVA.
% Outputs:
%   FC   - F-value for columns
%   dfc  - degree of freedom for columns
%
% Note: this function is inspired from rm_anova available at 
%       http://www.mathworks.se/matlabcentral/fileexchange/6874-two-way-rep
%       eated-measures-anova
%
% Example:
%   a = { rand(1,10) rand(1,10) rand(1,10) }
%   [FC dfc] = anova1rm_cell(a)
%   signifC = 1-fcdf(FC, dfc(1), dfc(2))
%
%   % for comparison 
%   [F1 F2 FI df1 df2 dfi] = anova1rm_cell(a);
%   F2
%
%   c = { rand(200,400,10) rand(200,400,10) };
%   [FC dfc] = anova2rm_cell(c) % computes 200x400 ANOVAs
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 2010

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

function [fA dfApair] = anova1rm_cell(data)

% compute all means and all std
% -----------------------------
a = length(data);
nd = myndims( data{1} );
sz = size( data{1} );
n  = size( data{1} ,nd);
AS = zeros([ sz(1:nd-1) a n ], 'single');
sq = zeros([ sz(1:nd-1) 1], 'single');

% only for paired stats
% ---------------------
for ind1 = 1:a
    switch nd
        case 1, AS(ind1,:)             = AS(ind1,:)             + reshape(data{ind1},[sz(1:nd-1) 1 n]);
        case 2, AS(:,ind1,:)           = AS(:,ind1,:)           + reshape(data{ind1},[sz(1:nd-1) 1 n]);
        case 3, AS(:,:,ind1,:)         = AS(:,:,ind1,:)         + reshape(data{ind1},[sz(1:nd-1) 1 n]);
        case 4, AS(:,:,:,ind1,:)       = AS(:,:,:,ind1,:)       + reshape(data{ind1},[sz(1:nd-1) 1 n]);
        case 5, AS(:,:,:,:,ind1,:)     = AS(:,:,:,:,ind1,:)     + reshape(data{ind1},[sz(1:nd-1) 1 n]);
        case 6, AS(:,:,:,:,:,ind1,:)   = AS(:,:,:,:,:,ind1,:)   + reshape(data{ind1},[sz(1:nd-1) 1 n]);
        case 7, AS(:,:,:,:,:,:,ind1,:) = AS(:,:,:,:,:,:,ind1,:) + reshape(data{ind1},[sz(1:nd-1) 1 n]);
        otherwise     error('Dimension not supported');
    end;
    sq = sq + sum(data{ind1}.^2,nd);
end;
dimA = nd+1;
dimB = nd;

A = sum(AS,dimA); % sum across columns, so result is 1xs row vector
S = sum(AS,dimB); % sum across columns, so result is 1xs row vector
T = sum(sum(S,dimB),dimA); % could sum either A or B or S, choice is arbitrary

% degrees of freedom
dfA = a-1;
dfAS = (a-1)*(n-1);

% bracket terms (expected value)
expA  = sum(A.^2,dimB)./n;
expS  = sum(S.^2,dimA)./a;
expAS = sum(sum(AS.^2,dimB),dimA);
expT  = T.^2 / (a*n);

% sums of squares
ssA   = expA - expT;
ssAS  = expAS - expA - expS + expT;

% mean squares
msA   = ssA / dfA;
msAS  = ssAS / dfAS;

% f statistic
fA = msA ./ msAS;
dfApair  = [dfA dfAS];

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
