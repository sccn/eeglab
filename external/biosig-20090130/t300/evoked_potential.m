function R = evoked_potential(fn,CHAN,t1,t2,EventTyp)
% EVOKED_POTENTIAL estimates evoked potentials (EP's)
%
%  R = EVOKED_POTENTIAL(filename, CHAN, t_start, t_end,EventTyp)
%  R = EVOKED_POTENTIAL(s, HDR, t_start, t_end,EventTyp)
%     filename  filename
%     CHAN      channel selection; default: 0 (all)
%     t_start   start time in seconds (relative to trigger time point)
%     t_end     end time in seconds relative to trigger 
%     EventTyp  [optional]
%
%  The trigger information must be available in the biosig file. 
%  The EP is calculated for each selected channel, if classlabels 
%  are available, the EP is calculated for each class
% 
%  The 

%	$Id: evoked_potential.m,v 1.1 2009-01-30 06:04:43 arno Exp $
%	Copyright (C) 2005,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 3 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.


%% 
%[S,HDR]=matload(fn,CHAN);
%S(S> 400)=NaN;S(S<-400)=NaN;
if ischar(fn)
	[S,HDR]=sload(fn,CHAN);
	if (CHAN==0)
		CHAN=1:HDR.NS;
	end; 	
	HDR.Classlabel = HDR.EVENT.TYP;
	HDR.TRIG = HDR.EVENT.POS;
elseif isnumeric(fn)
	S = fn; 
	HDR = CHAN;
end;		
%HDR = sopen(fn,'r',CHAN);[S,HDR]=sread(HDR);HDR = sclose(HDR); 
%S = diff(S); 
%[B,A]=butter(5,[450 900]/HDR.SampleRate*2);S = filtfilt(B,A,S);  
%[B,A]=fir1(20,[450 900]/HDR.SampleRate*2);S = filtfilt(B,A,S);  

if ~isfield(HDR,'TRIG') 
        error('trigger information is missing')
end;        
if ~isfield(HDR,'Classlabel');
        HDR.Classlabel = ones(size(HDR.TRIG)); 
end;        

if nargin<4,
        t1 = 2; t2 = 4;
end;

R0 = []; se=[];m=[];
CL = unique(HDR.Classlabel(:))';
if nargin>4,
	CL = EventTyp; 
end; 

t1 = floor(t1*HDR.SampleRate);
t2 = ceil(t2*HDR.SampleRate);
t  = t1:t2;
for cl = 1:length(CL), 
	if ~isfield(HDR,'Classlabel') || isempty(HDR.Classlabel)
		trig = HDR.EVENT.POS(HDR.EVENT.TYP==CL(cl)); 
	else	
		trig = HDR.TRIG(HDR.Classlabel==CL(cl));
	end;
	if 0, 
		sz = [length(CHAN),length(t),length(ix)];
		s  = repmat(NaN,sz);
		for k = 1:sz(3), 
			ix1 = ix(k)+t1;
			ix2 = ix(k)+t2;
			if (ix1>0) & (ix2<size(S,1)), 	
				%%%  FIXME: include also partial trials %%% 
				s(:,:,k) = S(ix1:ix2,:)';
			end;
		end;
	else 	 	
        	[s,sz] = trigg(S,trig,floor(t1*HDR.SampleRate),ceil(t2*HDR.SampleRate),15); 
        	N(cl)  = length(trig);
        	[se(:,:,cl), m(:,:,cl)] = sem(reshape(s,sz),3);
        end;  
        %RES = statistic(reshape(s,sz),3); 
        RES = statistic(center(s,2),3); 
        R0.SUM(:,:,cl) = RES.SUM';
        R0.N(:,:,cl)   = RES.N';
        R0.SSQ(:,:,cl) = RES.SSQ';
end; 
R0.datatype = 'MEAN+STD';
R0.T = t/HDR.SampleRate;
%R0.trigger = 3*HDR.SampleRate; 

R = R0; 
if all(size(CHAN)>1), 
elseif (CHAN==0)
	R.Label = HDR.Label;
	if isfield(HDR,'ELEC')
		R.ELEC = HDR.ELEC;
	end;
elseif all(CHAN>0)
	R.Label = HDR.Label(CHAN,:); 
	if isfield(HDR,'ELEC')
		R.ELEC.XYZ = HDR.ELEC.XYZ(CHAN,:);
	end; 	
end;

if nargout>0, return; end

R0.MEAN = R0.SUM ./ R0.N;			% mean 
R0.SSQ0	= R0.SSQ - real(R0.SUM).*real(R0.MEAN) - imag(R0.SUM).*imag(R0.MEAN);	% sum square of mean removed

sz = [size(R.SUM),1];
R.SUM = reshape(R.SUM,[sz(1),prod(sz(2:3))])'; 
R.N   = reshape(R.N,  [sz(1),prod(sz(2:3))])'; 
R.SSQ = reshape(R.SSQ,[sz(1),prod(sz(2:3))])'; 

R.MEAN 	= R.SUM./R.N;			% mean 
R.MSQ  	= R.SSQ./R.N;;			% mean square
R.RMS  	= sqrt(R.MSQ);			% root mean square
%R.SSQ0	= R.SSQ-R.SUM.*R.MEAN;		% sum square of mean removed
R.SSQ0	= R.SSQ - real(R.SUM).*real(R.MEAN) - imag(R.SUM).*imag(R.MEAN);	% sum square of mean removed
    n1 	= max(R.N-1,0);			% in case of n=0 and n=1, the (biased) variance, STD and SEM are INF
R.VAR  	= R.SSQ0./n1;	     		% variance (unbiased) 
R.STD  	= sqrt(R.VAR);		     	% standard deviation
R.SEM  	= sqrt(R.SSQ0./(R.N.*n1)); 	% standard error of the mean


if nargout>0, return; end

for k1 = 1:sz(2),
        for k2=1:sz(3),
                nf(k1,k2)=subplot(sz(2),sz(3),(k1-1)*sz(3)+k2); 
        end;
end;
if all(size(CHAN>1)), 
elseif (CHAN==0)
	R.Label = HDR.Label; 
elseif all(CHAN>0) 
	R.Label = HDR.Label(CHAN,:); 
end;
plota(R,nf)

for k=1:size(nf,1),
	set(nf(k),'YLabel',HDR.Label(CHAN(k),:))
end;	
for k=1:size(nf,2),
%	set(nf(k),'Title',sprintf('Class #))
end;	


