% circle() - draw a circle in the current Matlab axes
%
% Usage:
%   >> [linehandles fillhandle] = circle(X,Y,radius,colorfill,coloredge,oriangle,...
%                                         endangle,dashed,thickness,segments);
% Inputs:
%   X, Y       - circle center coordinates 
%   radius     - circle radius. Can be a vector of 2 values, one
%                for the x dimension and one the y dimension.
%   colorfill  - circle fill color (default:0=none)
%   coloredge  - circle edge color (default:black; 0:no edge)  
%   oriangle   - starting angle (default:0 degrees)
%   endangle   - ending angle (default:360 degrees)
%   dashed     - 0=no, 1=yes (default: no)
%   thickness  - thickness of boarder (default: line default)
%   segments   - number of line segments (default:50)
%
% Outputs:
%   linehandles - handle to the lines of the circle object
%   fillhandle  - handle to the interior of the circle object

% Author: Arnaud Delorme, CNL / Salk Institute, 2001

% This program is free software; you can redistribute it and/or modify it.  
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

function [h, h2] = circle( X, Y, radius, colorfill, coloredge, oriangle, endangle, dashed, thickness, segments);

if nargin < 3
	error('Not enough arguments. Type help circle');
	return;
end;
if nargin < 4
	colorfill = 0;
end;
if nargin < 5
	coloredge = 'k';		
end;
if nargin < 6
	oriangle = 0;		
end;
if nargin < 7
	endangle = 360;		
end;
if nargin < 8
	dashed = 0;
end;
if nargin < 9
	thickness = 0;
end;
if nargin < 10
	segments = 50;
end;
if any(radius <= 0)
    return;
end;

A = linspace(oriangle/180*pi, endangle/180*pi, segments);

% draw surface
% ------------
if any(colorfill)
	A = linspace(oriangle/180*pi, endangle/180*pi, segments);
	h2 = patch( X + cos(A)*radius(1), Y + sin(A)*radius(end), zeros(1,segments), colorfill);
	set(h2, 'FaceColor', colorfill); 
	set(h2, 'EdgeColor', 'none'); 
end;		

% draw lines
% ----------
if dashed
	compt=0;	
	for i=1:3:segments-2
		compt = compt+1;
		h(compt) = line( X + cos(A(i:i+1))*radius(1), Y + sin(A(i:i+1))*radius(end));
	end;
else
	h = line( X + cos(A)*radius(1), Y + sin(A)*radius(end));
end;
set(h, 'Color', coloredge); 

if thickness ~= 0
	set(h, 'Linewidth', thickness); 
end;
return;
