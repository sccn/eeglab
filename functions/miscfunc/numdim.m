% numdim() - estimate a lower bound on the (minimum) number of discrete sources 
%            in the data via their second-order statistics.
% Usage:
%   >> num = numdim( data );
%
% Inputs:
%   data   - 2-D data (nchannel x npoints)
%
% Outputs:
%   num    - number of sources (estimated from second order measures)
%
% References:
%   WACKERMANN, J. 1996. Beyond mapping: estimating complexity 
%   of multichannel EEG recordings. Acta Neurobiologiae 
%   Experimentalis, 56, 197-208.
%   WACKERMANN, J. 1999. Towards a quantitative characterization 
%   of functional states of the brain: from non-linear methodology 
%   to the global linear description. International Journal of 
%   Psychophysiology, 34, 65-80.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 23 January 2003

% Copyright (C) 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function lambda = numdim( a )
    
    if nargin < 1
        help numdim;
        return;
    end;
    
% Akaike, Identification toolbox (linear identification)

    a = a';
    b = a'*a/100; % correlation
    [v d] = eig(b);
    %det(d-b); % checking
    
    l = diag(d);
    l = l/sum(l);
    lambda = real(exp(-sum(l.*log(l))));
    
    return;
    
   
    % testing by duplicating columns
    a = rand(100,5)*2-1;
    a = [a a];
    numdim( a )
