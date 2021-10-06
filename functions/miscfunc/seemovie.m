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
	ntimes = -10;    % default to playing forward|backward endlessly
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
