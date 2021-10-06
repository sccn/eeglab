% chancenter() - recenter cartesian X,Y,Z channel coordinates
%
% Usage:  >> [x y z newcenter] = chancenter(x,y,z,center); 
%
% Optional inputs:
%    x,y,z     = 3D coordinates of the channels
%    center    = [X Y Z] known center different from [0 0 0]
%                [] will optimize the center location according
%                to the best sphere. Default is [0 0 0].
%
% Note: 4th input gui is obsolete. Use pop_chancenter instead.
%
% Authors: Arnaud Delorme, Luca Finelli & Scott Makeig SCCN/INC/UCSD,
%          La Jolla, 11/1999-03/2002 
%
% See also: spherror(), cart2topo()

% Copyright (C) 11/1999 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 3-16-00 improved help message -sm
% 1-25-02 put spherror subfunction inside chancenter -ad
% 1-25-02 help pops-up if no arguments -ad
% 01-25-02 reformated help & license -ad 
% 02-13-02 now center fitting works, need spherror outside chancenter -lf
% 02-14-02 radii are squeezed of squeeze in to fit topoplot circle -lf
% 03-31-02 center fitting is optional
% 04-01-02 automatic squeeze calculation -ad & sm
 
function [ x, y, z, newcenter, optim] = chancenter( x, y, z, center, gui)

optim = 0;

if nargin<4
    help chancenter
    return;
end

if nargin > 4 && gui
    error('Chancenter: 4th input'' gui'' is obsolete. Use pop_chancenter instead');
else 
	if isempty(center)
		optim = 1;
		center = [0 0 0];
	end
end

options = {'MaxFunEvals',1000*length(center)};
x = x - center(1);  % center the data
y = y - center(2);
z = z - center(3);
radius = (sqrt(x.^2+y.^2+z.^2));   % assume xyz values are on a sphere
if ~isempty(radius)
     wobble = std(radius);              % test if xyz values are on a sphere
else wobble = [];
end
fprintf('Radius values: %g (mean) +/- %g (std)\n',mean(radius),wobble);
newcenter = center;

if ~isempty(wobble) && wobble/mean(radius) > 0.01 && optim==1
	% Find center
	% ----------------------------------------------
	fprintf('Optimizing center position...\n');
	kk=0;
	while wobble/mean(radius) > 0.01 && kk<5
        try
    		newcenter = fminsearch('spherror',center,struct(options{:}),x,y,z);
        catch
            try
        		newcenter = fminsearch('spherror',center,options,x,y,z);
            catch
            	newcenter = fminsearch('spherror',center,[], [], x,y,z);
            end
        end
		nx = x - newcenter(1);  % re-center the data
		ny = y - newcenter(2);
		nz = z - newcenter(3);
		nradius = (sqrt(nx.^2+ny.^2+nz.^2));   % assume xyz values are on a sphere
		newobble = std(nradius);   
		if newobble<wobble
			center=newcenter;
			fprintf('Wobble too strong (%3.2g%%)! Re-centering data on (%g,%g,%g)\n',...
					100*wobble/mean(radius),newcenter(1),newcenter(2),newcenter(3))
			x = nx;  % re-center the data
			y = ny;
			z = nz;
			radius=nradius;
			wobble=newobble;
			kk=kk+1;
		else
			newcenter = center;
			kk=5;
		end
	end
	fprintf('Wobble (%3.2g%%) after centering data on (%g,%g,%g)\n',...
			100*wobble/mean(radius),center(1),center(2), ...
			center(3))
	%else
	%  fprintf('Wobble (%3.2g%%) after centering data on (%g,%g,%g)\n',...
	%              100*wobble/mean(radius),center(1),center(2),center(3))
end
