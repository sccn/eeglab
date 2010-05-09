% kurt() - return kurtosis of input data distribution
%
% Usage:
%   >> k=kurt(data)
%
% Algorithm:
%   Calculates kurtosis or normalized 4th moment of an input data vector
%   Given a matrix, returns a row vector giving the kurtosis' of the columns
%   (Ref: "Numerical Recipes," p. 612)
%
% Author: Martin Mckeown, CNL / Salk Institute, La Jolla, 10/2/96

% Copyright (C) Martin Mckeown, CNL / Salk Institute, La Jolla, 7/1996
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

% 2/28/97 - made to return separate kurtosis estimates of columns -Scott Makeig
% 01-25-02 reformated help & license, added links -ad 

function [k] = kurt(data)

[r,c]=size(data);
if r==1,
	kdata = data';  % if a row vector, make into a column vector
    r = c;
else
    kdata = data;
end
%fprintf('size of kdata = [%d,%d]\n',size(kdata,1),size(kdata,2));

mn = mean(kdata);              % find the column means
diff = kdata-ones(r,1)*mn;     % remove the column means
dsq = diff.*diff;              % square the data

k =  (sum(dsq.*dsq)./std(kdata).^4)./r - 3;

