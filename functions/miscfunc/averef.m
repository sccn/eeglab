% averef() - convert common-reference EEG data to average reference
%
% Usage:
%   >> data = averef(data);
%   >> [data W] = averef(data, W);
%   >> [data W S] = averef(data, W, S);
%
% Inputs & outputs:
%   data - 2D data (chans,frames*epochs) matrices of EEG or MEG data
%   W    - ICA weigth matrix
%   S    - ICA sphere matrix
%
% Note: the weight martix returned by ICA also have to be average
%       referenced:
%       because ICAACT = W*DATA, DATA = INV(W)*ICAACT 
%       as R=averef matrix, average reference data is R*DATA
%       so R*DATA = (R*INV(W))*ICA and NEW_W = INV(R*INV(W))
%
% Authors: Scott Makeig and Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 1999 

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 1999 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% 12/16/99 Corrected denomiator on the suggestion of Ian Nimmo-Smith, Cambridge UK
% 01-25-02 reformated help & license -ad 

function [data, W, S] = averef(data, W, S)

if nargin<1
  help averef
  return
end
chans = size(data,1);
if chans < 2 
  help averef
  return
end

% avematrix = eye(chans)-ones(chans)*1/chans;
% data = avematrix*data; % implement as a matrix multiply
% else (faster?)

data = data - ones(chans,1)*sum(data)/chans;

% treat optional ica parameters
if nargin == 2
	winv = pinv(W);
	avematrix = eye(chans)-ones(chans)*1/chans;
	W = pinv(avematrix*winv);
end;
if nargin == 3
	winv = pinv(W*S);
	avematrix = eye(chans)-ones(chans)*1/chans;
	W = pinv(avematrix*winv);
	S = eye(size(W,1), size(W,1));
end;
