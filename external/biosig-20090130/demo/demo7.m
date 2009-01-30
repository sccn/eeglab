function y = demo7(arg1)
% Demo7 - MULTIVARIATE AUTOREGRESSIVE Analysis: 
%
%      1) Simulates a MVAR process
%      2) Estimatess MVAR parameters
%      3) Displays the PDC for the original parameters and its estimates
%
%   demo7(k)
%       k       if k is a skalar between 1 to 5: use simulation k from Baccala et al. [1] 
%   demo7(eeg)
%       demonstrates the simulation of Kus et al. [3]        
%
% see also: BACCALA2001, MVAR, MVFILTER, PLOTA, demo/demo7
%
% Reference(s):
%  [1] Baccala LA, Sameshima K. (2001)
%       Partial directed coherence: a new concept in neural structure determination.
%       Biol Cybern. 2001 Jun;84(6):463-74. 
%  [2] M. Kaminski, M. Ding, W. Truccolo, S.L. Bressler, Evaluating causal realations in neural systems:
%	Granger causality, directed transfer functions and statistical assessment of significance.
%	Biol. Cybern., 85,145-157 (2001)
%  [3] R. Kus, M. Kaminski, K.J.Blinowska, Determination of EEG Activity Propagation - 
%       Pairwise vs. Multichannel Estimate. IEEE Trans. Biomedical
%       Engineering 51(9) 1501-1510 (Sep 2004);
%  [4] http://biosig.sf.net/

%	$Id: demo7.m,v 1.1 2009-01-30 06:04:39 arno Exp $
%	Copyright (C) 1999-2005 by Alois Schloegl <a.schloegl@ieee.org>
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the  License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% MULTIVARIATE AUTOREGRESSIVE Analysis: 
if (nargin<1) | all(size(arg1)==1),
        k = arg1;
        [a{1},a{2},a{3},a{4},a{5}] = baccala2001;       
        AR0 = a{k};        % select parameter set
        M = size(AR0,1);
        x = randn(M,10000);
else
        % Kus et al 2004
        M = 7; 
        eeg = arg1(:)'; 
        x = rms(eeg)/2*randn(7,length(eeg)); 
        
        x(1,:) = x(1,:)+eeg;
        AR0 = 0.8*[0,1,0,1,0,0,0; 0,0,1,0,0,0,0; 0,0,0,0,0,0,0; 0,0,0,0,1,1,0;0,0,0,0,0,0,0; 0,0,0,0,0,0,0; 0,0,0,0,0,0,0; ]';
        B = eye(M); 
end


% Simulated MVAR prosses can be produced with 
y = mvfilter(eye(M),[eye(M),-AR0],x);  

% Estimate AR parameters
[AR,RC,PE] = mvar(y',3);

% The PDF and the DTF can be displayed with the following functions
X.A = [eye(M),-AR0]; X.B = eye(M); X.C = eye(size(X.A,1));
X.datatype = 'MVAR';
X.SampleRate = 100;
X0 = X;

ffields = {'GPDC','logS','phase','COH','iCOH','pCOH','DTF','PDC','PDCF','dDTF','ffDTF','GGC','iSpectrum','rSpectrum','Af'};

%ffields = {'phase','dT'};
for k=1:length(ffields);
	figure(k)
	plota(X0,ffields{k});
%	figure(2)
%	plota(X1,ffields{k});
	drawnow;
	pause

end;
return;

figure(1)
plota(X,'PDC')         
if exist('suptitle','file')
        suptitle('PDC of true MVAR parameters');
end;
% The PDF and the DTF can be displayed with the following functions
X.A = [eye(M),-AR]; X.B = eye(M); X.C  = PE(:,(1-M:0)+end);
X.datatype = 'MVAR';

figure(2)
plota(X,'DTF')         
if exist('suptitle','file')
        suptitle('DTF of MVAR estimates');
end;
return;
figure(3)
plota(X,'COH')         
if exist('suptitle','file')
        suptitle('COH of MVAR estimates');
end;

figure(4)
plota(X,'SPECTRUM')         
if exist('suptitle','file')
        suptitle('Auto- and Cross-Spectra from MVAR estimates');
end;






