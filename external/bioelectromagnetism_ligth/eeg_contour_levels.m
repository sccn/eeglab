function [Clevels] = eeg_contour_levels( step, data )

% eeg_contour_levels - Defines contour levels relative to zero.
%
% Useage: [Clevels] = eeg_contour_levels( step, data )
%
%           step     - the step size
%           data     - array of values
%           Clevels  - array of values at steps
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:51 $

% Licence:  GNU GPL, no implied or express warranties
% History:  Created 07/2001, Darren.Weber_at_radiology.ucsf.edu
%           - developed under IDL and recoded for matlab
%
% BUG:      - may need to generate one more
%             negative contour level because matlab doesn't
%             fill the min contour level.
%           - 04/2002, Darren.Weber
%             modified lines 28 to 45 to attempt a fix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Clevels = [];

%Ensure step is a positive value
step = abs(step);


% Calculate min step range of voltage values
minstep = 0;
if min(data) < 0,
    while minstep - step > min(data),
        minstep = minstep - step;
    end
else
    while minstep < min(data),
        minstep = minstep + step;
    end
end

% Calculate max step range of voltage values
maxstep = 0;
if max(data) > 0
    while maxstep + step < max(data),
        maxstep = maxstep + step;
    end
else
    while maxstep > max(data),
        maxstep = maxstep - step;
    end
end



% Check for unusual conditions
if (minstep > maxstep),
    fprintf('...eeg_contour_levels, minstep > maxstep, revise step.\n');
    if minstep > min(data) & minstep < max(data),
        Clevels = minstep;
    end
    if maxstep > min(data) & maxstep < max(data),
        Clevels = maxstep;
    end
    return
end
if and(maxstep == 0, minstep == 0),
    fprintf('...eeg_contour_levels, step too large for min/max of data.\n');
    return
end



% Define Clevels according to min/max step relative to zero

if and(minstep < 0.0, maxstep > 0.0)
    
	% Define number of negative levels, given step and range
	nlevels = round( ( 0 - minstep ) / step );
	% Define column (Narray) to hold (nlevels) from -step to minstep
    for a = 0:(nlevels -1), Narray(1,a+1) = (-1 * step) - (a * step); end

	% Define number of positive levels, given step and range
	plevels = round( maxstep / step );
	% Define column (Parray) to hold (plevels) from zero to maxstep
	for a = 0:plevels, Parray(1,a+1) = step * a; end

	Clevels = [fliplr(Narray), Parray];
    return
end

if and(minstep < 0.0, maxstep <= 0.0)

	% Define number of negative levels, given step and range
	nlevels = round( ( maxstep - minstep ) / step );
	% Define (Narray) to hold nlevels from maxstep to minstep
    % (max[data] may not be in Narray).
    for a = 0:nlevels, Narray(1,a+1) = maxstep - (step * a); end

	Clevels = fliplr(Narray);
    return
end

if and(minstep >= 0.0, maxstep > 0.0)

	% Define number of positive levels, given step and range
	plevels = round( ( maxstep - minstep ) / step );
	% Define (Parray) to hold (plevels) from minstep to maxstep,
    % (min[data] may not be in Parray).
	for a = 0:plevels, Parray(1,a+1) = minstep + (step * a); end

	Clevels = Parray;
    return
end
