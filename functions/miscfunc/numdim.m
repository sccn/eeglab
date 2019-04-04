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

function lambda = numdim( a )
    
    if nargin < 1
        help numdim;
        return;
    end
    
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
