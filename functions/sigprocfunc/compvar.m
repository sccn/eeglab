% compvar()   - project selected components and compute the variance of
%               the original data they account for.
%
% Usage:
%   >> [proj, variance] = compvar( data, wts_or_act, winv, components);
%
% Required Inputs:
% data        - 2-D (channels, points) or 3-D (channels, frames, trials)
%               data array.
% wts_or_act  - {sphere weights} cell array containing the ICA sphere 
%               and weights matrices. May also be a 2-D (channels, points) 
%               or 3-D (channels, frames, trials) array of component 
%               activations.
% winv        - inverse or pseudo-inverse of the product of the weights
%               and sphere matrices returned by the ICA decompnumsition,
%               i.e., inv(weights*sphere) or pinv(weights*sphere).
% components  - array of component indices to back-project
%
% Outputs:
%  proj       - summed back-projections of the specified components
%  pvaf       - percent variance of the data that the selected 
%               components account for (range: 100% to -Inf%).
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% $Log: not supported by cvs2svn $
% Revision 1.7  2007/08/17 16:49:55  scott
% reworked help msg; made 2nd output pvaf; rm'd mean of data
%
% Revision 1.6  2007/08/07 01:33:57  arno
% fix help message
%
% Revision 1.5  2003/10/09 22:35:49  arno
% correcting sphere*weights
%
% Revision 1.4  2003/01/09 18:29:34  arno
% debuging compvar
% absolute value instead of square
%
% Revision 1.3  2002/04/10 19:34:39  arno
% futher debuging
%
% Revision 1.2  2002/04/10 19:31:12  arno
% debugging (variable name error)
%
% Revision 1.1  2002/04/05 17:39:45  jorn
% Initial revision
%
% 01-25-02 reformated help & license -ad 
% 03-08-02 added the sphere and weight option -ad

function [ compproj, varegg ] = compvar( data, act, winv, compnums);

if nargin < 4
   help compvar;
   return;
end;   

data = reshape(data, size(data,1), size(data,2)*size(data,3));
squaredata  = sum(sum(data.^2));             % compute the grand sum-squared data

if iscell(act)
    sphere = act{1};
    weight = act{2};
    act = (weight(compnums,:)*sphere)*data;
end;

compproj   = winv(:,compnums)*act-data;      % difference between data and back-projection
squarecomp = sum(sum(compproj.^2));          % the summed-square difference
varegg     = 100*(1- squarecomp/squaredata); % compute pvaf of components in data
compproj   = compproj+data;                  % restore back-projected data

return;

