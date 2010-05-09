% superline() - dashed (moving) line between two points
%
% Usage:
%   >> [handles] = superline( X, Y, color, linewidth, offset, middle );
%
% Inputs:
%   X          - abscices of the points (2 abscices for the 2 points)
%   Y          - ordinates of the points (2 ordinates for the 2 points)
%   color      - color of the line (default black)
%   linewidth  - thickness of the line (same as line 'linewidth' property, default:1)   
%   offset     - change the offset of the dash (in between 0 and 1) (default:0)
%   middle     - put a small line at the middle
%
% Outputs:
%   handles    - handles of all the objects composing the line

% arno@salk.edu, Arnaud Delorme, CNL / Salk Institute, 2001

% This program is free software; you can redistribute it and/or
% modify it.  
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

function h = superline(X, Y, color, thickness, offset, middle);

if nargin < 2
	error('Not enough arguments. Type help superline');
	return;
end;
if nargin < 3
	color = 'k';
end;
if nargin < 4
	thickness = 1;
end;
if nargin < 5
	offset = 0;		
end;
if nargin < 6
	middle = 0;
end;

% ratio X/Y
% ---------
if abs(Y(1) - Y(2)) == 0
	ratioXY = 1;
else 
	ratioXY = abs(X(1) - X(2))/abs(Y(1) - Y(2));
end;
		
PADSTRIPE  = 0.8;
SIZESTRIPE = 0.2;
%PADSTRIPE  = 0.18;
%SIZESTRIPE = 0.09;

% original line vector
% --------------------
if ratioXY >= 0.05
	slope = (Y(2)-Y(1)) / (X(2)-X(1));
	slopenormalize = [1 slope]/sqrt(slope^2+1); % (1,slope) are the coordinate of the non unitary vector
else
	slopenormalize = sign(Y(2)-Y(1))*[0 1];
end;	
% reverse vector if necessary
if X(1) > X(2)
	slopenormalize = - slopenormalize;
end;	
vectpad    = slopenormalize*PADSTRIPE*abs(complex(X(1), Y(1)) - complex(X(2), Y(2)));
vectstripe = slopenormalize*SIZESTRIPE*abs(complex(X(1), Y(1)) - complex(X(2), Y(2)));
nblines    = 1/(PADSTRIPE+SIZESTRIPE);

% draw the discontinuous line
% ---------------------------
cont = 1;
compter = 1;
firstpass = 1;
posx1 = X(1);
posy1 = Y(1);
offsetshifted = mod( 0.5 + offset + SIZESTRIPE/2,1);

while cont == 1
	if firstpass == 1 
		if offsetshifted/nblines < SIZESTRIPE
			posx1 = posx1 + vectstripe(1)*(offsetshifted/nblines)/SIZESTRIPE; 
			posy1 = posy1 + vectstripe(2)*(offsetshifted/nblines)/SIZESTRIPE; 
			posx2 = posx1 + vectpad(1);
			posy2 = posy1 + vectpad(2);
		else
			posx2 = posx1 + vectstripe(1)*(offsetshifted/nblines-SIZESTRIPE)/SIZESTRIPE; %consider the offset from 0 to 1
			posy2 = posy1 + vectstripe(2)*(offsetshifted/nblines-SIZESTRIPE)/SIZESTRIPE; %consider the offset from 0 to 1
		end;	
		firstpass = 0;
		%posy1
		%posy2
	else
		posx2 = posx1 + vectpad(1);
		posy2 = posy1 + vectpad(2);
		%posy1
		%posy2
	end;
	%fprintf(['Abscice1 ' num2str(X(1)) '  Abscice2 ' num2str(X(1)) '  equal ' num2str(X(1) == X(2)) '\n']);
	if ratioXY < 0.1
		if Y(1) > Y(2)
			if posy2 <= Y(2)
				posx2 = X(2);
				posy2 = Y(2);
				cont = 0;
			end;
		else
			if posy2 >= Y(2)
				posx2 = X(2);
				posy2 = Y(2);
				cont = 0;
			end;
		end;		
	else	
		if X(1) > X(2)
			if posx2 <= X(2)
				posx2 = X(2);
				posy2 = Y(2);
				cont = 0;
			end;
		else
			if posx2 >= X(2)
				posx2 = X(2);
				posy2 = Y(2);
				cont = 0;
			end;
		end;
	end;			
	h(compter) = line([posx1 posx2], [posy1 posy2]);
	set(h(compter), 'color', color);
	set(h(compter), 'linewidth', thickness);
	posx1 = posx2+vectstripe(1);
	posy1 = posy2+vectstripe(2);
	compter = compter +1;
	if ratioXY < 0.1
		if Y(1) > Y(2)
			if posy1 <= Y(2)
				cont = 0;
			end;
		else
			if posy1 >= Y(2)
				cont = 0;
			end;
		end;		
	else	
		if X(1) > X(2)
			if posx1 <= X(2)
				cont = 0;
			end;
		else
			if posx1 >= X(2)
				cont = 0;
			end;
		end;
	end;			
end;

%patch(patchccordinatesX, patchccordinatesY, [0 0 0 0], color);

if middle
	hold on; h = plot((X(1)+X(2))/2, (Y(1)+Y(2))/2, '*k', 'markersize', thickness/2);
end;

return

% tests
figure;
for i=9:31
	fprintf('******************************* %d \n', mod(i/10, 1));
	superline( [ 0+i 1+i ] , [ 10 20 ], 'b', 5, mod(i/10, 1),1);  
end;	




