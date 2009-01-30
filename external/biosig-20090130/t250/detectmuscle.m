function [INI,S,E] = detectmuscle(S, iter, Mode)
% Muscle detection with inverse filtering
% Artifacts are indicated with NaN's. 
%
% [INI,S,E] = detectmuscle(S [, iter [,1]])
% [INI,S,E] = detectmuscle(S [, Fs, 2])
%
% iter		number of iterations [default:1]
% INI.MU	mean of S
% INI.InvFilter	coefficients of inverse filter 
% S		outlier replaced by NaN
% E		innan(E) indicates muscle artifact


%	$Id: detectmuscle.m,v 1.1 2009-01-30 06:04:43 arno Exp $
%	Copyright (C) 2003,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/
%
%    BioSig is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BioSig is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with BioSig.  If not, see <http://www.gnu.org/licenses/>.


if nargin<2,
        iter=[];
end;
if isempty(iter) iter = 1; end; 

% Muscle Detection
if Mode==1, 
	%% inverse filter   
	TH = 5; 
	[se,INI.MU] = sem(S);
	if TH*se<abs(INI.MU),
		[S] = center(S);
	end;

	while iter>0,
        	[AR,RC,PE]= lattice(S',10);
		INI.InvFilter = ar2poly(AR);
        	E = zeros(size(S));
        	for k = 1:size(S,2),
                	E(:,k) = filter(INI.InvFilter(k,:),1,S(:,k));
	        end;
        	INI.V  = std(E);
	        INI.TH = INI.V * TH; 
        	iter   = iter-1;
        
	        for k = 1:size(S,2),
        	        S(E(:,k)>(INI.TH(k)) | E(:,k)<(-INI.TH(k)),k) = NaN;
	        end;
	end; 
	% the following part demonstrates a possible correction 
	for k = 1:size(S,2),
                E(:,k) = filter(INI.InvFilter(k,:),1,S(:,k));
	end;
	% isnan(E), returns Artifactselse 

else
	Fs = iter; 
	% Criterion for bad gradient:
	% Maximal allowed voltage step / sampling point: 100.00 µV
	% Mark as bad before event: 100.00 ms
	% Mark as bad after event: 100.00 ms
	ix1 = [abs(diff(S,[],1))>100;zeros(1,size(S,2))];
	
	% Criterion for bad max-min:
	% Maximal allowed absolute difference: 150.00 µV
	% Interval Length: 100.00 ms
	% Mark as bad before event: 100.00 ms
	% Mark as bad after event: 100.00 ms
	
	ix2 = abs(filter([-1;zeros(Fs/10-1,1);1],1,S))>150;
	ix2 = [ix2(Fs/20+1:end,:);zeros(Fs/20+1,size(S,2))]; 

	% Criterion for bad amplitude:
	% Minimal allowed amplitude: -75.00 µV
	% Maximal allowed amplitude: 75.00 µV
	% Mark as bad before event: 100.00 ms
	% Mark as bad after event: 100.00 ms
	ix3 = abs(S)>75;
	
	ix = filtfilt(ones(Fs/10,1),1,double(ix1|ix2|ix3))>0;
	S(ix) = NaN;	
	INI=[]; 

end; 

if nargout<3, return; end;

