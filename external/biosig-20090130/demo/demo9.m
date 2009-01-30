
% DEMO2 demonstrates the use of the data set III from the BCI competition 2003 for 
%   The demo shows the offline analysis for obtaining a classifier and 
%   uses a jack-knife method (leave-one-trial out) for validation. 
%   AAR parameters are extracted 
%
%
% References: 
% [1] A. Schlögl, B. Kemp, T. Penzel, D. Kunz, S.-L. Himanen, A. Värri, G. Dorffner, G. Pfurtscheller.
%       Quality Control of polysomnographic Sleep Data by Histogram and Entropy Analysis.
%       Clin. Neurophysiol. 1999, Dec; 110(12): 2165 - 2170.
% [2] Alois Schlögl (2000)
%       The electroencephalogram and the adaptive autoregressive model: theory and applications
%       Shaker Verlag, Aachen, Germany, (ISBN3-8265-7640-3). 
% [3] Schlögl A., Neuper C. Pfurtscheller G.
%       Estimating the mutual information of an EEG-based Brain-Computer-Interface
%       Biomedizinische Technik 47(1-2): 3-8, 2002.
% [4] A. Schlögl, C. Keinrath, R. Scherer, G. Pfurtscheller,
%       Information transfer of an EEG-based Bran-computer interface.
%       Proceedings of the 1st International IEEE EMBS Conference on Neural Engineering, Capri, Italy, Mar 20-22, 2003. 
% [5] A. Schlögl, J. Kronegg, J.E. Huggins, S. G. Mason.
%       Evaluation criteria in BCI research.
%       (Eds.) G. Dornhege, J.R. Millan, T. Hinterberger, D.J. McFarland, K.-R.Müller,
%       Towards Brain-Computer Interfacing. MIT Press, p.327-342, 2007.
% [6] A. Schlögl, F.Y. Lee, H. Bischof, G. Pfurtscheller
%   	Characterization of Four-Class Motor Imagery EEG Data for the BCI-Competition 2005.
%   	Journal of neural engineering 2 (2005) 4, S. L14-L22
% [7] A. Schlögl, C. Brunner, R. Scherer, A. Glatz;
%   	BioSig - an open source software library for BCI research.
%   	(Eds.) G. Dornhege, J.R. Millan, T. Hinterberger, D.J. McFarland, K.-R. Müller;
%   	Towards Brain-Computer Interfacing, MIT Press, 2007, p.347-358. 
% [8] A. Schlögl, C. Brunner
%   	BioSig: A Free and Open Source Software Library for BCI Research.
%	Computer (2008, In Press)	

%	$Id: demo9.m,v 1.1 2009-01-30 06:04:39 arno Exp $
%	Copyright (C) 1999-2003,2006,2007,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/


% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'B0101T.gdf';
	
% Step 1a: load data ============================%
fprintf(1,'Step 1: Prepare data.\n',filename);
fprintf(1,'\ta: Load data file %s.\n',filename);
[s,HDR] = sload(filename);

CHAN = sort([strmatch('ECG',HDR.Label);strmatch('EKG',HDR.Label);strmatch('ecg',HDR.Label);strmatch('Ecg',HDR.Label)]);
if length(CHAN)<1,
	fprintf(stderr,'no ECG channel found\n');
	return; 
elseif length(CHAN)>1,	 
	fprintf(stderr,'more than one ECG channel is currently not supported.\n');
	return;
else 	
	s = s(:,CHAN); 	 
end; 	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Preprocessing and artifact processing ====================%
%   2a: QRS detection
% QRS-Detection
H2  = qrsdetect(s,HDR.SampleRate);
RRI = diff(H2.EVENT.POS)/H2.EVENT.SampleRate; 

%   2b: outlier detection and rejection
%RRI(RRI>2)==NaN; 
%RRI(RRI<0.5)==NaN; 

%   2c: scale data 
s = log(RRI); 

%    3b: Adaptive Autoregressive Parameters
fprintf(1,'\tb: Adaptive Autoregressive parameters (AAR).\n');
f2 = [];
MODE.MOP  = [1,12,0];
MODE.UC   = 0.0085; 
for ch    = 1:length(eegchan),
	X = tvaar(s, MODE.MOP, MODE.UC); 	% get initial values
	X = tvaar(s, X);		% AAR estimation
end;
f0 = exp(-X.AAR(:,1)./(1-sum(X.AAR(:,2:end),2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 5: Show results  ==================== %

plota(X, f0, '3D'); 
