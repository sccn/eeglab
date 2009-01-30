function [X] = criteria4asyncbci(D0, TRIG, ONTIME ,Fs)
% CRITERIA4ASYNCBCI provides an evaluation criterion of asychronous BCI. 
% based on the discussion at the BCI2005 meeting in Rensellearville.  
%
% X = CRITERIA4ASYNCBCI(D, TRIG, ONTIME, Fs)
% X = CRITERIA4ASYNCBCI(D, STATE, [], Fs)
%   D           detector output
%               can be continuous; D should contain as many traces (columns) states. 
%               [ max(STATES) < size(D,2) ]
%               if all values are lower than the Threshold=0, the 
%               no-control (NC) state is assumed. 
%               If D is discrete, it must contain only one column, but no AUC can 
%               estimated. 
%   TRIG        list of trigger events in seconds 
%               if TRIG is a cell-array of a list of trigger times, 
%                       TRIG{n} contains the list of trigger times for
%                       state n. The number of detector traces (columns)
%                       must be 1+number of states.
%   ONTIME      on-time [in seconds] of event-driven discrete control
%   Fs          sampling rate (default = 1Hz)
%   STATE       "true" state 
%
%   X.H         confusion matrix       
%   X.AUC       area-under-the-curve.
% 
% X = CRITERIA4ASYNCBCI(D, TRIG, ONTIME, Fs)
%               TRIG and ONTIME are used to construct the STATE
%
%   X.TPR       True Positive Ratio
%   X.FPR       False Positive Ratio / False alarm rate
%
% X = CRITERIA4ASYNCBCI(D, STATE, [], Fs)
%   D           State of detector output (discrete values, D=0 indicates no-control state) 
%   STATE       Target state 
% 
%
%
% 
%

% References: 
%


%    $Id: criteria4asyncbci.m,v 1.1 2009-01-30 06:04:50 arno Exp $
%    Copyright (C) 2005 by Alois Schloegl <a.schloegl@ieee.org>	
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


%%%%% INTERNAL VARIABLES %%% 
%% D0   detector output, continuous in time and magnitude, one trace for each state
%% D    classified transducer output 
%%              [tmp,D] = max(D0,[],2); 
%%              D(isnan(tmp)) = Nan;    % not classified 
%%              D(tmp<Threshold) = 0;   % NC state
%% STATE    target class; continuos in time, discrete states, same size than D
%%              STATE can be constructed from TRIG and ONTIME from
%%              event-driven discrete control. 
%%  
%%
%%

Mode = 1;       % basic level, NC and one STATE state
X.datatype = 'Criteria_Asychronous_BCI'; 

if nargin<4,
        Fs = 1; 
end;

ONTIME = round(ONTIME*Fs);

M1 = size(D0,2);         
if iscell(TRIG(1))      % multiple states
        M2 = length(TRIG); 
        if length(TRIG)~=size(D0,2),
                warning('number of states does not fit number of detector traces');     
        end;
        T = []; CL = [];
        for k = 1:length(TRIG)
                t  = TRIG{k}(:);
                T  = [T;t]; 
                CL = [CL; repmat(k,length(t),1)];
        end;  
        TRIG = T; 
        Mode = 2; 
else
        M2 = 1; 
        CL = ones(numel(TRIG),1);
end;

FLAG.CONTINOUS = M1>1; 
if ~FLAG.CONTINOUS 
        tmp = double(D0(~isnan(D0))); 
        FLAG.CONTINOUS  = ~all(round(tmp)==tmp);
end;
% classifies transducer output / applies threshold

if size(D0,2)>1,
        [tmp,D] = max(D0,[],2);
        D(isnan(tmp)) = NaN; 
        D(tmp<0)= 0;                     % D=0,1,2, ... indicate output states N, A, B, ... respectively,
elseif FLAG.CONTINOUS,
        D = D0>0; 
else 
        D = D0; 
        M1 = max(D);
end; 

if isempty(ONTIME )
        STATE = TRIG;  
else         
        [TRIG,ix] = sort(TRIG*Fs); 
        CL = CL(ix);    % classlabels 
        
        if any(diff(TRIG) < ONTIME )
                warning('overlapping detection window - ONTIME reduced');
                ONTIME = min(diff(TRIG))-1;
        end;
        %### OPEN QUESTION(s): 
        %###    is there a reasonable way to deal with overlapping windows ? 
        
        STATE  = zeros(size(D));
        d1 = real(D>0);
        for k = 1:length(TRIG)                
                d1(TRIG(k):TRIG(k)+ONTIME ) = NaN;       % de-select (mark with NaN) all samples within window
                STATE(TRIG(k):TRIG(k)+ONTIME )  = CL(k);     % generate target state from TRIG and ONTIME 
        end;
end;

        % signal detection theory applied on a sample-basis for multiple classes
[X.kappa, X.kapSD, H, z] = kappa(D(:), STATE(:));
X.AS.H   = H; 					% confusion matrix 


X.AS.AUC = NaN; 	
if FLAG.CONTINOUS %any(size(D0,2)==[1,M])
	for k = 1:M1,
    		X.AS.AUC(k) = auc(D0(:,k),STATE==k); 
	end;
end;

%%% still missing 


if any([M1,M2]>1), return; end; 
    % following code not ready for more than 1 control state


N0     = length(TRIG);                      % number of trigger events 
% Detection of True Positives
[d,sx] = trigg(D,TRIG,0,ONTIME );           % intervals where detection is counted as hit
d      = reshape(d,sx(2:3));               % and bring it in proper shape
TP     = sum(any(d>0,1));                  % true positives/hits 

N3     = sum(diff(D>0)>0);                  % total number of detections (JH) 


% Detection of False Positives
[FP1,N1] = sumskipnan(d1);              % false positive ratio on a sample-basis
FP3  = sum(diff(d1)>1);                 % number of detections

N2   = ceil((1+sum(~~diff(isnan(d))))/2);          % number of intervals between trigger events (with non-overlapping windows) 
d(isnan(d)) = [];
tmp  = diff(d)>0;
FP2  = sum(tmp);                       % false positives
FPR2 = FP2/N2;                         % false positive ratio 


%  suggestion by Steve. # is this correct? 
X.SM.TP  = TP; 
X.SM.TPR = TP/N0; 
X.SM.FP  = FP1; 
X.SM.FPR = FP1/N1;
X.SM.N   = N1; 

% Jane Huggins' suggestion 
X.JH.TP  = TP; 
X.JH.TPR = TP/N0; 
X.JH.FP  = FP3; 
X.JH.FPR = FP3/N3; 
X.JH.N   = N3; 
X.JH.HFdiff = max(0,X.JH.TPR-X.JH.FPR);        % ???

%### OPEN QUESTION: use of FPR1 or FPR2
X.TPR = TP/N0;                % true positive ratio
X.FPR = FP1/N1;                % or 
%X.FPR = FP2/N2;        % dismissed 
%X.FPR = FP3/N3; 
X.FPR = '?';





return;
% summary statistics
X.HFdiff = max(0,X.TPR-X.FPR);         % 
X.Dprime = norminv(X.TPR,0,1)-norminv(X.FPR,0,1);
tmp = X.TPR-X.FPR;
X.Aprime = .5+sign(tmp)*(tmp*tmp+abs(tmp))/(4*max(X.TPR,X.FPR)-4*X.TPR*X.FPR);
tmp2 = [X.TPR*(1-X.TPR),X.FPR*(1-X.FPR)];
X.Bsecond = sign(tmp)*diff(tmp2)/sum(tmp2);
