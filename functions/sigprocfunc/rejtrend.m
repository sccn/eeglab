% rejtrend() - detect linear trends in EEG activity and reject the  
%                  epoched trials based on the accuracy of the linear
%                  fit.
% Usage:
%   >> [rej rejE] = rejtrend( signal, winsize, maxslope, minR, step);
%
% Inputs:
%   signal     - 3 dimensional signal (channels x frames x trials)
%   winsize    - integer determining the number of consecutive points
%                for the detection of linear patterns
%   maxslope   - maximum acceptable absolute slope of the linear trend. If the slope
%                of the line fitted to a data epoch is greater than or equal to 
%                maxslope, that epoch is rejected (assuming a sufficient R^2 value, 
%                see minR below).
%   minR       - minimal R^2 (coefficient of determination between
%                0 and 1). The R^2 value reflects how well the data epoch is
%                approximated by a line. An epoch is not rejected unless its R^2
%                value is greater than minR.
%   step       - step for the window. Default is 1 point (2 points 
%                will divide by two the computation time)  
%
% Outputs:
%   rej        - rejected trials. Array with 0 or 1 for each trial.
%   rejE       - rejected rows of the rejected trials
%
% Algorithm:
%   Looked for all possible windows of size 'winsize' of each trial if 
%   the linear fit have minimum slope and R^2
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [rej, rejE] = rejtrend( signal, pointrange, maxslope, minstd, step);

% This is to avoid divide-by-zero and machine errors.
SST_TOLERANCE = 1000*pointrange*1.1921e-07;

if nargin < 3
	help rejtrend
	return;
end;	
if nargin < 5
	step = pointrange;
end;	
[chans pnts trials] = size(signal);
rejE   = zeros( chans, trials);

x = linspace( 1/pointrange, 1, pointrange );
%waitbarhandle = waitbar(0,'rejtrend.m Please wait...');
for c = 1:chans
	for t = 1:trials 
		for w = 1:step:(pnts-pointrange+1)
			y = signal(c, [w:w+pointrange-1], t);
		   	coef = polyfit(x,y,1);   		
			if abs(coef(1)) >= maxslope
			   	ypred = polyval(coef,x);   % predictions
			   	dev = y - mean(y);          % deviations - measure of spread
			   	SST = sum(dev.^2);          % total variation to be accounted for
                if SST < SST_TOLERANCE      % make sure SST is not too close to zero
                    SST = SST_TOLERANCE;
                end
			   	resid = y - ypred;              % residuals - measure of mismatch
			   	SSE = sum(resid.^2);           % variation NOT accounted for
                Rsq = 1 - SSE/SST;             % percent of error explained
				if Rsq > minstd
					rejE( c, t ) = 1;
                end
				% see the page 	http://www.facstaff.bucknell.edu/maneval/help211/fitting.html
            end
        end
    end
    %waitbar(c/chans);
end
%close(waitbarhandle);			
rej = max( rejE, [], 1);

return;
