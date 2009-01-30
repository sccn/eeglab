function [K1,K2,K3] = criteria4momentarybci(TO,Fs,trig,ERW)
% CRITERIA4MOMENTARYBCI evaluates the output of momentary self-paced BCI  
%
% [K1,K2,K3] = criteria4momentarybci(TO,Fs,TRIG,ERW) 
%
% Input: 
% 	TO transducer output: 0 is no control state, i>0 is i-th control state
% 	Fs sampleing rate
% 	TRIG 	trigger information - its a vector in case of a single IC state
%		othewise TRIG{K} are the triggers for the K-th IC state
% 	ERW 	expected response window in seconds, default=[-0.5,+0.5]
%
% output:
% 	K1 	result the rule-based (heuristic) approach [1] 
% 	K2 	results of the sample-by-sample approach
%	K3	result by Fatourechi's Method 
%
% 	K1.HF  		hf-difference 
% 	K{123}.kappa 	Cohen's kappa coefficient
% 	K{123}.H 	Confusion matrix 
% 
% References: 
% [1] Jane E. Huggins, Michael T. McCann, …
% 	Comparison of Evaluation Metrics for an ECoG-based Momentary BCI.
% 	https://ctools.umich.edu/access/content/attachment/5fa4908e-eaba-4e67-8eee-d92169c70726/kappa_vs_HF_windowed2_jeh3.doc
% [2] M Fatourechi, R K Ward, and G E Birch
%	A self-paced brain–computer interface system with a low false positive rate
%	J. Neural Eng. 5: 9–23, 2008.

%    $Id: criteria4momentarybci.m,v 1.1 2009-01-30 06:04:50 arno Exp $
%    Copyright (C) 2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    This is part of the BIOSIG-toolbox http://biosig.sf.net/
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


if nargin<4
	ERW = [-.5,.5];
end; 	

EUI = zeros(size(TO));
if iscell(trig)
	CL = [];
	T = [] 
	for K = 1:length(trig)
		T = [T; trig{K}(:)];
		CL = [CL; repmat(K,length(trig{K}),1)];
		for k=1:length(TRIG)
			if any(EUI((TRIG(k)+ERW(1)*Fs) : (TRIG(k)+ERW(2)*Fs)))
				warning('overlapping ERW segments');
			end; 	
			EUI((TRIG(k)+ERW(1)*Fs) : (TRIG(k)+ERW(2)*Fs))=K; 
		end;
	end; 
	[T,ix]=sort(T); 
	CL = CL(ix);
else 
	TRIG = trig; 
	CL   = ones(size(TRIG));
	for k= 1:length(TRIG)
		if any(EUI((TRIG(k)+ERW(1)*Fs) : (TRIG(k)+ERW(2)*Fs)))
			warning('overlapping ERW segments');
		end; 	
		EUI((TRIG(k)+ERW(1)*Fs) : (TRIG(k)+ERW(2)*Fs)) = CL(k); 
	end; 
end; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Jane's rules for the evaluation 

EV  = zeros(length(TO),2);  % default value is TN 

for k=1:length(TRIG)

% 1) The parameters ERWstart and ERWend are chosen to define the ERW 
%	around each momentary event label. For the present comparison, 
%	ERWstart is .5 seconds before the trigger, and ERWend is .5 seconds after the trigger.

	ix = max(1,TRIG(k)+ERW(1)*Fs):min(TRIG(k)+ERW(2)*Fs,length(TO)); 

% 2) The first detection within the range of ERWstart to ERWend 
%	relative to the trigger is considered a TP.

	F = find(TO(ix));
	if ~isempty(F)
		EV(ix(F(1)),1)=CL(k);  	
		EV(ix(F(1)),2)=TO(ix(F(1)));  	
		
% 3) Any additional detection(s) within the ERW after the first is considered a FP.

		EV(ix(F(2:end)),1) = 0; 		
		EV(ix(F(2:end)),2) = TO(ix(F(2:end)));  		

% 4) If there are no detections within the window, a single 
%	sample in the center of the window is marked FN.

	else 
		EV(TRIG(k),:)=[1,0]; 
	end; 

end;

% 5) Any detection outside of an ERW is considered a FP.

	ix = find(~EUI);
	F  = find(TO(ix));
	EV(ix(F),1) = 0;  % FP 		
	EV(ix(F),2) = TO(ix(F));  % FP 		
		
% 6) All other points (both inside and outside of ERWs) are considered TN.

	% initial value of EV	


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% summary statistics: compute kappa and hf-difference
K1 = kappa(EV(:,1),EV(:,2)); 

if length(K1.H)==2,
	%% HF difference = TP/(TP+FN) - FP/(TP+FP)
	K1.HF = K1.H(2,2)/sum(K1.H(2,:)) - K1.H(1,2)/sum(K1.H(:,2));  
end; 	


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sample-based evaluation 
%% generate padding window - make padding window the same size than ECW
T1 = TO;   
i0 = find(TO);
for k = i0(:)',
	cl = TO(k); 
	ix = max(1,k + ERW(1)*Fs) : min(k + ERW(2)*Fs,length(TO));
	if any((T1(ix)~=cl) & (T1(ix)>0))
		warning('overlapping padding segments with different IC states');
	end; 
	T1(ix) = cl; 
end; 
K2 = kappa(EUI,T1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mehrdad's method on "TNs based on switch output activation" [2]

EV  = [zeros(length(TO),1),TO(:)];  % default value is TN or FP  	 
for k = 1:length(TRIG)
	ix = max(1,TRIG(k)+ERW(1)*Fs):min(TRIG(k)+ERW(2)*Fs,length(TO));

	% overwrite values within ECW 
	if any(TO(ix))	% TP
		EV(ix(1),1)=1;  	
		EV(ix(1),2)=1;  	
		EV(ix(2:end),:)=NaN;	% NaN's are ignored, therefore just 1 count per ECW
	else 		% FN
		EV(ix(1),1)=1;  	
		EV(ix(1),2)=0;  	
		EV(ix(2:end),:)=NaN;	% NaN's are ignored, therefore just 1 count per ECW
	end; 	
end;
EV = EV(~any(isnan(EV),2),:);		
K3 = kappa(EV(:,1),EV(:,2)); 

