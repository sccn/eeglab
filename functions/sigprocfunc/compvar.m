% compvar() - project selected components and compute the variance of
%             the original signal they account for.
%
% Usage:
%   >> [proj, variance] = compvar( rawdata, ica_act, winv, components);
%
% Inputs:
%   rawdata    - 2D(channels x points) or 3D(channels x frames x trials)
%              data array.
%   ica_act    - 2D(channels x points) or 3D(channels x frames x trials)
%              ica activity array. It can also be a cell array with the
%              sphere and the weight matrix { sphere weight } that will
%              be used to recompute the ICA activity.
%   frames     - number of points per epoch.
%   winv       - inverse of (sphere * weights) matrices returned by the
%              ica function. It represent the component distribution of
%              activity across the electrodes.
%   components - array of components to project.
%
% Outputs:
%   proj       - projection of the components
%   variance   - variance of the original signal that the selected 
%              components account for.
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

function [ compproj, varegg ] = compvar( sig, act, winv, compos);

if nargin < 4
   help compvar;
   return;
end;   

sig = reshape(sig, size(sig,1), size(sig,2)*size(sig,3));
squaresig  = sum(sum(sig.^2));

if ~iscell(act)
    compproj   = winv(:,compos)*act(compos,:)-sig;
else
    weight = act{2};
    sphere = act{1};
    acttmp = (weight(compos,:)*sphere)*sig;
    compproj   = winv(:,compos)*acttmp-sig;
end;

squarecomp = sum(sum(compproj.^2));
varegg     = 1- squarecomp/squaresig;
compproj   = compproj+sig;

%meaneeg( s ) = mean(sum(compproj.*compproj))/(size(compproj,1)-1);

return;

% old version
% -----------
for s = compos
    % compute projection
    % ------------------
    fprintf('%d ',s);         % construct single-component data matrix

    %compproj = (winv(:,s)*act(s,:))./sig(:,:);
    %[I] = find(compproj > 0);
    %meaneeg( s ) = mean(mean(compproj(I)));

    compproj = winv(:,s)*act(s,:)-sig(:,:);
	squarecomp = sum(compproj.*compproj);
 	squaresig  = sum(sig(:,:).*sig(:,:));
    varegg( s ) = squarecomp / squaresig;

    %meaneeg( s ) = mean(sum(compproj.*compproj))/(size(compproj,1)-1);
end;
fprintf('\n');
