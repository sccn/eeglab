% rejtrend() - detect linear trends in EEG activity and reject the  
%                  epoched trials based on the accuracy of the linear
%                  fit.
% Usage:
%   >> [rej rejE] = rejtrend( signal, winsize, minslope, minR, step);
%
% Inputs:
%   signal     - 3 dimensional signal (channels x frames x trials)
%   winsize    - integer determining the number of consecutive points
%                for the detection of linear patterns
%   minslope   - minimal absolute slope of the linear trend of the 
%                activity for rejection
%   minR       - minimal R^2 (coefficient of determination between
%                0 and 1)
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
% 01-25-02 reformated help & license -ad 

function [rej, rejE] = rejtrend( signal, pointrange, minslope, minstd, step);

if nargin < 3
	help rejtrend;
	return;
end;	
if nargin < 5
	step = 1;
end;	
[chans pnts trials] = size(signal);
rejE   = zeros( chans, trials);

% normalize each row input
% ------------------------
%signal = signal(:,:);
%signal = (signal-mean(signal,2)*ones(1,size(signal,2)))./ (std(signal,0,2)*ones(1,size(signal,2)));
%signal = reshape(signal, chans, pnts, trials);

%fprintf('finding stable low variability regions of %d consecutive time points...\n', pointrange);	
x = linspace( 1/pointrange, 1, pointrange );
for c = 1:chans
	%fprintf('%d ', c);
	for t = 1:trials 
		for w = 1:step:(pnts-pointrange+1)
			y = signal(c, [w:w+pointrange-1], t);
		   	coef = polyfit(x,y,1);   		
			if abs(coef(1)) > minslope
			   	ypred = polyval(coef,x);   % predictions
			   	dev = y - mean(y);          % deviations - measure of spread
			   	SST = sum(dev.^2);          % total variation to be accounted for
			   	resid = y - ypred;          % residuals - measure of mismatch
			   	SSE = sum(resid.^2);        % variation NOT accounted for
			   	Rsq = 1 - SSE/SST;          % percent of error explained
				if Rsq > minstd
					rejE( c, t ) = 1;
				end;
				% see the page 	http://www.facstaff.bucknell.edu/maneval/help211/fitting.html
			end;			
			%rejE( c, t ) = rejE( c, t ) | all( abs(detrend(signal(c, [w:w+pointrange-1], t))) < minstd);
			%rejE( c, t ) = rejE( c, t ) | all( abs(signal(c, [w:w+pointrange-1], t)) < minstd);
		end;
	end;
end;
%fprintf('\n');					
rej = max( rejE, [], 1);

return;
