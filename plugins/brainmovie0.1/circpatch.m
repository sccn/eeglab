% circpatch() - draw a dashed arc between two points
%
% Usage:
%   >> [handles] = circpatch( X, Y, curv, color, linewidth, segments, offset );
%
% Inputs:
%   X          - abscissae of the points (2 abscices for the 2 points)
%   Y          - ordinates of the points (2 ordinates for the 2 points)
%   curv       - curvature (0 = no curvature, 1 = round, -1 = round in 
%                the other direction)
%   color      - line color (default black)
%   linewidth  - line thickness (same as line 'linewidth' property, default:1)   
%   segments   - number of segments composing each continuous part of the arc (default:50)
%   offset     - change the offset of the dash (range: [0,1]) (default:0)
%   middle     - put a small line at the middle of the arc
%
% Outputs:
%   handles    - handles of all the objects composing the arc

% arno@salk.edu, Arnaud Delorme, CNL / Salk Institute, 2001

% This program is free software; you can redistribute it and/or modify it.  
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

function h = circpatch( X, Y, circfactor, color, thickness, segments, offset, middle );

PADSTRIPE  = 0.8; % percentage
SIZESTRIPE = 0.2; % percentage

if nargin < 3
	error('Not enough arguments. Type help circpatch');
	return;
end;
if nargin < 4
	color = 'k';
end;
if nargin < 5
	thickness = 1;
end;
if nargin < 6
	segments = 50;
end;
if nargin < 7
	offset = 0;		
end;
if nargin < 8
	middle = 0;		
end;

% middle of X and Y
% -----------------
middleXY = [(X(1)+X(2))/2 (Y(1)+Y(2))/2];

% normal to this line
% -------------------
diffXY = [ (X(2)-X(1))/2 (Y(2)-Y(1))/2];
normdiffXY = [diffXY(2) -diffXY(1)];

% third coordinate
% ----------------
fisrtPoint  = [ X(1) Y(1) ];
secondPoint = [ X(2) Y(2) ];
thirdPoint  = middleXY + normdiffXY*circfactor;

% compute the two lines equation 
% ------------------------------
a1 = normdiffXY(2)/normdiffXY(1);
b1 = thirdPoint(2)-thirdPoint(1)*a1;

middle13 = [(fisrtPoint(1)+thirdPoint(1))/2 (fisrtPoint(2)+thirdPoint(2))/2];
a2 = [1 (fisrtPoint(2)-thirdPoint(2))/(fisrtPoint(1)-thirdPoint(1))];
a2 = [a2(2) -a2(1)]; % take the normal
a2 = a2(2)/a2(1);
b2 = middle13(2)-middle13(1)*a2;

% solve the equation (the center of the circle is the solution)
% -------------------------------------------------------------
centreX = (b1-b2)/(a2-a1);
centreY = a1*centreX + b1;
centreY = a2*centreX + b2;
radius  = sqrt((fisrtPoint(1)-centreX)^2+(fisrtPoint(2)-centreY)^2);

%circle( fisrtPoint(1), fisrtPoint(2), 0.05);
%circle( secondPoint(1), secondPoint(2), 0.05);
%circle( thirdPoint(1), thirdPoint(2), 0.05);

% compute the angle of the points
% -------------------------------
angle1 = atan( (fisrtPoint(2)-centreY)/(fisrtPoint(1)-centreX) )/pi*180;
angle2 = atan( (secondPoint(2)-centreY)/(secondPoint(1)-centreX) )/pi*180;
if ((fisrtPoint(1)-centreX) < 0) % negative cosinus
	angle1 = 180+angle1;
end;
if ((secondPoint(1)-centreX) < 0)
	angle2 = angle2+180;
end;
if angle1 < 0 angle1 = angle1+360; end;
if angle2 < 0 angle2 = angle2+360; end;

if circfactor<0 % which side
	tmp = angle1;
	angle1 = angle2;
	angle2 = tmp;
end;	
if angle2 < angle1 angle2 = angle2+360; end;

% if offset exist, drawstripes
% ----------------------------  
if offset ~= -1
	offsetshifted = mod(0.5 + offset + SIZESTRIPE/2,1);
	vectpad    = PADSTRIPE*abs(angle2-angle1);
	vectstripe = SIZESTRIPE*abs(angle2-angle1);
	nblines    = 1/(PADSTRIPE+SIZESTRIPE);
	currentangle = angle1;
	cont =1;
	h = [];
	firstpass = 1;
	while cont
		if firstpass == 1
			if offsetshifted/nblines < SIZESTRIPE
				currentangle = currentangle + vectstripe*(offsetshifted/nblines)/SIZESTRIPE; 
				targetangle = currentangle + vectpad;
			else
				targetangle = currentangle + vectstripe*(offsetshifted/nblines-SIZESTRIPE)/SIZESTRIPE; %consider the offset from 0 to 1
			end;	
			firstpass = 0;
		else
			targetangle  = currentangle + vectpad;
		end;
		if targetangle >= angle2
			targetangle = angle2;
			cont = 0;
		end;
		hh = circle( centreX, centreY, radius, 0, color, currentangle, targetangle, 0, thickness, ceil((targetangle-currentangle)/4*radius) );
		h = [h hh];

		currentangle = targetangle+vectstripe;
		if currentangle >= angle2
			cont = 0;
		end;	
	end;	
else
	h = circle( centreX, centreY, radius, 0, color, angle1, angle2, 0, thickness, segments );
end;

if middle
	x = centreX + cos((angle1+angle2)/2/180*pi)*radius;
	y = centreY + sin((angle1+angle2)/2/180*pi)*radius;
	hold on; h = plot( x, y, '*k', 'markersize', thickness/2);
end;

return;

clf;
for i=9:31
	fprintf('******************************* %d \n', mod(i/10, 1));
	circpatch( [ 0+i 1+i ] , [ 10 20 ], 0.5, 'b', 5, 10, mod(i/10, 1),1);  
end;	
