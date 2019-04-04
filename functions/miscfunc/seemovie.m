% seemovie() - see an EEG movie produced by eegmovie()
%
% Usage: >> seemovie(Movie,ntimes,Colormap)
%
% Inputs:
%         Movie    = Movie matrix returned by eegmovie()
%         ntimes   = Number of times to display {0 -> -10}
%                    If ntimes < 0, movie will play forward|backward
%         Colormap = Color map returned by eegmovie() {0 -> default}
%
% Author: Scott Makeig & Colin Humphries, CNL / Salk Institute, 6/3/97
%
% See also: eegmovie()

% Copyright (C) 6/3/97 Scott Makeig & Colin Humphries, CNL / Salk Institute, La Jolla CA
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

% 10-31-97 changed test for user-defined colormap -ch & sm
% 1-8-98 added '\n' at end, improved help msg -sm
% 01-25-02 reformated help, added link -ad 

function seemovie(Movie,ntimes,Colormap)

fps = 10;   % projetion speed (requested frames per second)

if nargin<1
   help seemovie
   return
end
if nargin<3
    Colormap = 0;
end
if nargin<2
	ntimes = -10;    % default to playing foward|backward endlessly
end
if ntimes == 0
	ntimes = -10;
end

clf
axes('Position',[0 0 1 1]);
if size(Colormap,2) == 3 % if colormap user-defined
	colormap(Colormap)
else
	colormap([jet(64);0 0 0]);    % set up the default topoplot color map
end

if ntimes > 0,
 fprintf('Movie will play slowly once, then %d times faster.\n',ntimes);
else 
 fprintf('Movie will play slowly once, then %d times faster forwards and backwards.\n',-ntimes);
end
if abs(ntimes) > 7
	fprintf('   Close figure to abort:  ');
end

%%%%%%%%%%%%%%%%%%%%%%%% Show the movie %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

movie(Movie,ntimes,fps);  

%%%%%%%%%%%%%%%%%%%%%%%%    The End     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if abs(ntimes) > 7
	fprintf('\n');
end
