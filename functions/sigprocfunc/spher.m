% spher() - return the sphering matrix for given input data
%
% Usage:
%        >> sphere_matrix = spher(data);
%
% Reference: T. Bell (1996)

% Copyright (C) S. Makeig CNL / Salk Institute, La Jolla CA 7-17-97
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

function sphere = spher(data)


if nargin<1 || size(data,1)<1 
  help spher
  return
end

sphere = 2.0*inv(sqrtm(cov(data'))); % return the "sphering" matrix
