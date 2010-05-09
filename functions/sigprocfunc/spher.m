% spher() - return the sphering matrix for given input data
%
% Usage:
%
%        >> sphere_matrix = spher(data);
%
% Reference: T. Bell (1996) - - -
%

% S. Makeig CNL / Salk Institute, La Jolla CA 7-17-97

function sphere = spher(data)


if nargin<1 | size(data,1)<1 
  help spher
  return
end

sphere = 2.0*inv(sqrtm(cov(data'))); % return the "sphering" matrix
