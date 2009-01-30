function [S,HDR,N] = denoise(D,ny,nx)
%  DENOISE yields the regression coefficients for 
%  correcting EOG artifacts in EEG recordings.
%  
%  The correction of a single record is obtained like this:      
%   [S] = DENOISE(S)
%   [S] = DENOISE(FileName)
%  removes bias, and 50/60Hz.      
% 
%   [S] = DENOISE(..., SigChanIdx)
%  removes bias, and 50/60Hz on the Signal channels (indicated by vector
%  SigChanIdx) only. 
%
%   [S] = DENOISE(..., ListSigChan, ListNoiseChan)
%  uses the Noise channels for correction, too. 
% 
%  Input parameters:
%    filename   signal data file. 
%    S          signal data, each column represents a channel
%    SigChanIdx list of signal channels; 0 [default] indicates all channels
%    NoiseChanIdx  list of noise channels use for correction 
%  Output parameters
%    S          corrected signal data
%       
% see also: 
%
% Reference(s):
%  R.T.Behrens, L.L. Scharf, Signal Processing Applications of Oblique
%  Projection Operators. IEEE Trans Signal Processing, 42(6), 1994,
%  1413-1424. 
%

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

%	$Id: denoise.m,v 1.1 2009-01-30 06:04:43 arno Exp $
%	(C) 2006 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/


if nargin<2,
        ny=0; 
end; 
if nargin<3,
        nx = [];
end; 
if ischar(D),
        [D,HDR] = sload(D);
else
        HDR.NS = size(D,2);
        HDR.SampleRate = 1000; 
end

if (ny==0), ny = 1:size(D,2); end;      % default: correct all channels. 

D = D(1:1000,:);

a = speye(HDR.NS+1);
if size(nx,1)==1,  % list of noise channel 
        nx  = nx(:);
        
else    % noise channels are defined through rereferencing (e.g. bipoloar channels)
        tmp = HDR.NS+(1:size(nx,2)); 
        a(2:size(nx,1)+1,tmp+1) = nx;
        nx = tmp';
end;

[m,n]=size(D); 
is_missing = any(isnan(D),2);
N  = [D(:,nx),ones(size(D,1),1)];   % noise components, bias
if any(is_missing )      % if there are any missing samples
        N(find(is_missing ),:) = 0;
        N = [N,is_missing ];
end;
if isfield(HDR,'SampleRate'); 
        NotchFreq = 50;
        N = [N,real(exp((i*2*pi*NotchFreq/HDR.SampleRate)*[1:m]')*[1,i])];
        %N = [N,real(exp((i*2*pi*NotchFreq/HDR.SampleRate)*[1:m]'))*[1,-1]];
else    
        warning('Line interference 50/60 Hz not removed');
end;
E = eye(m) - N*pinv(N);   % P_S|

H = D(:,ny); 
H(is_missing,:) = 0;   % get rid of NaN's
%E_HS = H*((H'*E*H)\H'*E);  % E_HS
E_HS = H*pinv(E*H);

S = E_HS*[D(:,ny),N];     % use of D ensures that NaN's remain in the data. !!! What about NaN's in the noise? 
%N = E_HS*N;

