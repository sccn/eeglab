% mapcorr() - Find matching rows in two matrices and their corrs.
%             Uses the Hungarian (default), VAM, or maxcorr assignment methods.
%             (Follow with matperm() to permute and sign x -> y).
%
%             Finds correlation of maximum common subset of channels (using
%             channel location files to match channel labels.) Thus, number
%             of channels can differ in x and y.
%
% Usage: 
%   >> [corr,indx,indy,corrs] = matcorr(x,y,ch1,ch2);
%   >> [corr,indx,indy,corrs] = matcorr(x,y,ch1,ch2,rmmean,method,weighting);
%
% Inputs:
%   x     = first input matrix. Row are difference components and columns
%           the channels for these components.
%   y     = matrix with same number of columns (channels) as x
%   ch1   = channel locations file for x
%   ch2   = channel locations file for y
% 
% Optional inputs:
%   rmmean    = When present and non-zero, remove row means prior to correlation 
%               {default: 0}
%   method    = Method used to find assignments.
%               0= Hungarian Method - maximize sum of abs corrs {default: 2}
%               1= Vogel's Assignment Method -find pairs in order of max contrast 
%               2= Max Abs Corr Method - find pairs in order of max abs corr 
%               Note that the methods 0 and 1 require matrices to be square.
%   weighting = An optional weighting matrix size(weighting) = size(corrs) that 
%               weights the corrs matrix before pair assignment {def: 0/[]->ones()}
% Outputs:
%   corr  = a column vector of correlation coefficients between 
%           best-correlating rows of matrice x and y
%   indx  = a column vector containing the index of the maximum 
%           abs-correlated x row in descending order of abs corr 
%           (no duplications)
%   indy  = a column vector containing the index of the maximum 
%           abs-correlated row of y in descending order of abs corr 
%           (no duplications)
%   corrs = an optional square matrix of row-correlation coefficients
%           between matrices x and y
%
% Note: outputs are sorted by abs(corr)
%
% Authors: Scott Makeig & Sigurd Enghoff, SCCN/INC/UCSD, La Jolla, 11-30-96 

% Copyright (C) 11-30-96 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 2007/03/20 02:33:48  arno, Andreas fix
% 2003/09/04 23:21:06  scott, changed default matching method to Max Abs Corr
% 04-22-99 Re-written using VAM by Sigurd Enghoff, CNL/Salk
% 04-30-99 Added revision of algorthm loop by SE -sm
% 05-25-99 Added Hungarian method assignment by SE
% 06-15-99 Maximum correlation method reinstated by SE
% 08-02-99 Made order of outpus match help msg -sm
% 02-16-00 Fixed order of corr output under VAM added method explanations, 
%          and returned corr signs in abs max method -sm
% 01-25-02 reformated help & license, added links -ad 

% Uses function hungarian.m

function [corr,indx,indy,corrs] = mapcorr(x,y,ch1,ch2,rmmean,method,weighting)
%
if nargin < 4
  help matcorr
  return
end

if nargin < 6
	method = 2; % default: Max Abs Corr - select successive best abs(corr) pairs
end

[m,n] = size(x);
[p,q] = size(y);
m = min(m,p);

if m~=n | p~=q
   if nargin>5 & method~=2
     fprintf('matcorr(): Matrices are not square: using max abs corr method (2).\n');
   end
   method = 2; % Can accept non-square matrices
end 


if nargin < 5 | isempty(rmmean)
  rmmean = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if rmmean
  x = x - mean(x')'*ones(1,n); % optionally remove means
  y = y - mean(y')'*ones(1,n);
end

for i = 1:m
    for j = 1:p
        corrs(i,j) = compcorr(x(i,:)',ch1,y(j,:)',ch2);
    end
end

%dx = sum(x'.^2);
%dy = sum(y'.^2);
%dx(find(dx==0)) = 1;
%dy(find(dy==0)) = 1;
%corrs = x*y'./sqrt(dx'*dy);



if nargin > 6 & ~isempty(weighting) & norm(weighting) > 0,
  if any(size(corrs) ~= size(weighting))
    fprintf('matcorr(): weighting matrix size must match that of corrs\n.')
    return
  else
	corrs = corrs.*weighting;
  end
end

cc = abs(corrs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch method
case 0
	ass = hungarian(-cc);   % Performs Hungarian algorithm matching

	idx1 = sub2ind(size(cc),ass,1:m);
	[dummy idx2] = sort(-cc(idx1));
	corr = corrs(idx1);
	corr = corr(idx2)';
	indy = [1:m]';
	indx = ass(idx2)';
	indy = indy(idx2);

case 1                      % Implements the VAM assignment method
	indx = zeros(m,1);
	indy = zeros(m,1);
	corr = zeros(m,1);

	for i=1:m,
		[sx ix] = sort(cc);  % Looks for maximum salience along a row/column
		[sy iy] = sort(cc'); % rather than maximum correlation.
		[sxx ixx] = max(sx(end,:)-sx(end-1,:));
		[syy iyy] = max(sy(end,:)-sy(end-1,:));

		if sxx == syy
			if sxx == 0 & syy == 0
        	   [sxx ixx] = max((sx(end,:)-sx(end-1,:)) .* sx(end,:));
    	       [syy iyy] = max((sy(end,:)-sy(end-1,:)) .* sy(end,:));
			else
				sxx = sx(end,ixx); % takes care of identical vectors
				syy = sy(end,iyy); % and zero vectors
			end
		end

		if sxx > syy
			indx(i) = ix(end,ixx);
			indy(i) = ixx;
		else
			indx(i) = iyy;
			indy(i) = iy(end,iyy);
		end
		cc(indx(i),:) = -1;
		cc(:,indy(i)) = -1;
	end

	i = sub2ind(size(corrs),indx,indy);
	corr = corrs(i);

	[tmp j] = sort(-abs(corr)); % re-sort by abs(correlation)
	corr = corr(j);
	indx = indx(j);
	indy = indy(j);

case 2                       % match successive max(abs(corr)) pairs
	indx = zeros(size(cc,1),1);
	indy = zeros(size(cc,1),1);
	corr = zeros(size(cc,1),1);

	for i = 1:size(cc,1)
		[tmp j] = max(cc(:));
		% [corr(i) j] = max(cc(:));
		[indx(i) indy(i)] = ind2sub(size(cc),j);
        corr(i) = corrs(indx(i),indy(i));
		cc(indx(i),:) = -1; % remove from contention
		cc(:,indy(i)) = -1;
	end

otherwise
	error('Unknown method');
end


function corr = compcorr(a1,ch1,a2,ch2)

n1 = length(a1);
n2 = length(a2);

%if n1 < n2
    cnt = 0;
    corr = 0;
    for i = 1:n1
        for j = 1:n2
            if strcmp(ch1(i).labels,ch2(j).labels)
                cnt = cnt+1;
                b1(cnt,1) = a1(i);
                b2(cnt,1) = a2(j);
            end
        end
    end 
 b1 = b1 / norm(b1);
 b2 = b2 / norm(b2);
 corr = b1'*b2;
