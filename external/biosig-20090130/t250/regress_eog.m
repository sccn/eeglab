function [R,S] = regress_eog(D,ny,nx)
%  REGRESS_EOG yields the regression coefficients for 
%  correcting EOG artifacts in EEG recordings.
%  
%  The correction of a single record is obtained like this:      
%   [R,S2] = regress_eog(S1, EL, OL)
%   [R,S2] = regress_eog(filename, EL, OL)
%   [R,S2] = regress_eog(filename)
%	OL = IDENTIFY_EOG_CHANNELS(filename)
%	EL are all remaining channels
%       
%  Corrected data is obtained through
%   [R] = regress_eog(covm(S1,'E'), EL, OL)
%   S2 = S1 * R.r0;    % without offset correction
%   S2 = [ones(size(S1,1),1),S1] * R.r1;    % with offset correction
%  
%  S1   recorded data
%  EL   list of eeg channels: those channels will be corrected   
%  OL   eog/ecg channels. 
%       if OL is a vector, it represents the list of noise channels 
%       if OL is a matrix, OL derives the noise channels through rereferencing. 
%          This is useful if the EOG is recorded monopolar, but the bipolar EOG 
%          should be used for for artefact reduction (because the global EEG should remain), 
%          One can define OL = sparse([23,24,25,26],[1,1,2,2],[1,-1,1,-1]) 
%	   resulting in two noise channels defined as bipolar channels #23-#24 and #25-#26
%  R.r1, R.r0    rereferencing matrix for correcting artifacts with and without offset correction
%  R.b0	coefficients of EOG influencing EEG channels
%  S2   corrected EEG-signal      
%
% see also: IDENTIFY_EOG_CHANNELS, SLOAD
%
% Reference(s):
% [1] Schlogl A, Keinrath C, Zimmermann D, Scherer R, Leeb R, Pfurtscheller G. 
%	A fully automated correction method of EOG artifacts in EEG recordings.
%	Clin Neurophysiol. 2007 Jan;118(1):98-104. Epub 2006 Nov 7.
% 	http://dx.doi.org/10.1016/j.clinph.2006.09.003
%       http://www.dpmi.tugraz.at/~schloegl/publications/schloegl2007eog.pdf

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
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

%	$Id: regress_eog.m,v 1.1 2009-01-30 06:04:43 arno Exp $
%	(C) 1997-2007,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

if ischar(D),
        [D,H] = sload(D);
        C = covm(D,'E');
        if (nargin<3)
        	nx = identify_eog_channels(H); 
       	end; 
        if (nargin<2)
               	ny = find(~any(nx,2)); 
       	end; 
        
elseif size(D,1)==size(D,2)
        C = D; 
        H.NS = size(C,2)-1;

else
        H.NS = size(D,2);
        C = covm(D,'E');
end

R.datatype = 'ArtifactCorrection_Regression';
R.signalchannel = ny;
R.noise_channel = nx;
R.Mode = 1; % correction of the mean 
R.Mode = 0; % correction without mean 

a = speye(H.NS+1);
if size(nx,1)==1,  % list of noise channel 
        nx  = nx(:);
else    % noise channels are defined through rereferencing (e.g. bipoloar channels)
        tmp = H.NS+(1:size(nx,2)); 
%	nx(isnan(nx)) = 0;
        a(2:size(nx,1)+1,tmp+1) = nx;
        if rank(full(nx)) < size(nx,2),
                fprintf(2,'Warning REGRESS_EOG: referencing matrix is singular! \n');
        end; 
        C  = a'*C*a;
        nx = tmp';
end;
ny = ny(:);

r0 = speye(H.NS);
r1 = sparse(2:H.NS+1,1:H.NS,1);

b0 = -inv(C([1; nx+1], [1; nx+1])) * C([1; nx+1], ny+1);
r0(nx,ny) = b0(2:end,:);       % rearrange channel order
r1([1;1+nx], ny) = b0;       % rearrange channel order

R.b0 = b0;	% correction coefficients, useful for visualization
R.r0 = a(2:end,2:end)*r0;
R.r1 = a*r1;

if size(D,1)==size(D,2),
        % R = R;         
        
elseif (nargout>1)
        % S = D * R.r0;   % without offset correction 
        S = [ones(size(D,1),1),D] * R.r1; % with offset correction
end

