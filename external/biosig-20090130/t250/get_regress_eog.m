function [h0, s00] = get_regress_eog(fn,Mode)
% GET_REGRESS_EOG tries to obtain the regression coefficients
%    for EOG correction. According to [1], some extra recordings 
%    with large eye movements (i.e. EOG artifacts) are needed. 
%    GET_REGRESS_EOG tries to identify this data. 
% 
%   hdr = get_regress_eog(file)
%   hdr = get_regress_eog(file, Mode)
% 
% INPUT: 
%   file	filename which should be corrected.
%		usually, the eye movements are stored in a different file.
%    		Some lab-specific heuristics is used to identify the file with the eye movements.
%   Mode	'REG'	[default] regression with one or two bipolar EOG channels [1]
%		'REG+CAR' regression and common average reference
%			removes 2 bipolar + averaged monopolar EOG
%		'REG+PCA' regression and PCA, 
%			removes 3 "EOG" components
%		'REG+ICA' regression + ICA [3]
%			removes 3 "EOG" components 
%		'PCA-k'	removes the k-largest PCA components, k must be a positive integer
%		'ICA-k'	removes the k-largest ICA components [3], k must be a positive integer
%       	others like 'NGCA-k','TDSEP-k','TDSEP3-k','TDSEP1','FFDIAG'
%		'msec'	same as PCA-3, modified (without averaging) MSEC method [2]
%		'bf-'	beamformer, assume zero-activity reference electrode
%		'bf+'	beamformer, take into account activity of reference electrode 
%		'Hurst'	ICA components selected by the method of[5]
%		'Joyce2004' [6]
%		'Barbati2004' [7]	
%		'Meinecke2002' [8]	
%
%	The following modifiers can be combined with any of the above	
%		'FILT###-###Hz'  filtering between ### and ### Hz. ### must be numeric
%		'Fs=###Hz'  downsampling to ### Hz, ### must be numeric 
%		'x'  	2nd player of season2 data
%
% OUTPUT:
%   hdr.REGRESS.r0 	correction coefficients
%
%   The EOG correction will be applied to the channels CHAN with any of these commands: 
%   	HDR = sopen(file,'r',hdr.REGRESS.r0(:,CHAN)); [s,HDR]=sread(HDR); HDR=sclose(HDR);
%   	[s,HDR] = sload(file,hdr.REGRESS.r0(:,CHAN)); 
%   	[s,HDR] = sload(file,CHAN,'EOG_CORRECTION','ON'); 
%
% See also: SLOAD, IDENTIFY_EOG_CHANNELS, BV2BIOSIG_EVENTS, REGRESS_EOG 
%
% Reference(s):
% [1] Schlogl A, Keinrath C, Zimmermann D, Scherer R, Leeb R, Pfurtscheller G. 
%	A fully automated correction method of EOG artifacts in EEG recordings.
%	Clin Neurophysiol. 2007 Jan;118(1):98-104. Epub 2006 Nov 7.
% 	http://dx.doi.org/10.1016/j.clinph.2006.09.003
%       http://www.dpmi.tugraz.at/~schloegl/publications/schloegl2007eog.pdf
% [2] Berg P, Scherg M.
%	A multiple source approach to the correction of eye artifacts.
%	Electroencephalogr Clin Neurophysiol. 1994 Mar;90(3):229-41.
% [3] JADE algorithm, Jean-François Cardoso.
% [4] Boudet S., Peyrodie L., P Gallois, C Vasseur, 
% 	Filtering by optimal projectsion and application to automatic artifact removal from EEG
% 	Signal Processing 87 (2007) 1987-1992. 
% [5] Vorobyov and Cichocki (2002)
%	Blind noise reduction for multisensory signals using ICA and subspace 
%	filtering, with application to EEG analysis.
%	Biol Cybern. 2002 Apr;86(4):293-303.
% [6] C.A. Joyce, I.F. Gorodnitsky, M.Kutas
%	Automated removal of eye movement and blink artifats from EEG data using blind component separation.
%	Psychobiology, 41 (2004), 313-325
% [7] Barbati et al (2004) 
% [8] Frank Meinecke, Andreas Ziehe, Motoaki Kawanabe, and Klaus-Robert Müller.
%	A Resampling Approach to Estimate the Stability of One-Dimensional or Multidimensional Independent Components.
%	IEEE Transactions on Biomedical Engineering, 49(12):1514-1525, 2002.
% [9] Blanchard G., Kawanabe M., Sugiyama M., Spokoiny V., Muller K.-R. (2006). 
%	In search of non-gaussian components of a high-dimensional distribution. 
%	Journal of Machine Learning Research 7, 247-282.
% [10] Kawanabe M., Sugiyama M., Blanchard G, Müller K.-R. (2007) 
%	A new algorithm of non-Gaussian component analysis with radial kernel functions
%	Annals of the Institute of Statistical Mathematics, 59(1):2007
% [11] K.H. Ting, P.C.W. Fung, C.Q.Chang, F.H.Z.Chan
%	automatec correction of artifact from single-trial event-related potentials bz blind separation  using second order statistics only.
%	Medical Engineering & Physics, 28, 780-794 (2006)

%	$Id: get_regress_eog.m,v 1.1 2009-01-30 06:04:43 arno Exp $
%	Copyright (C) 2006,2007 by Alois Schloegl 
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 2 of the License, or (at your option) any later version.

	% extract information from BBCI data for EOG correction
	if nargin<3,
		CHAN = 0 ; 
	end; 	
	if nargin<2,
		Mode = ''; 
	end; 	

	[p,f,e] = fileparts(fn); 

	%%%=== search for eye-movement recordings (in order to obtain clean EOG artifacts) ===%%%
	%	an heuristic approach is to support lab-specific standards
	% 	currently, the convention of the GrazBCI and the BBCI Lab are supported
	
	f0 = dir(fullfile(p,['arte*',e])); % BBCI files
	LAB_ID = 'BBCI'; 
	if length(f0)>1,
		% hack to ignore "arte*eog_mono" data
		ix = regexp(strvcat({f0.name}),'eog_mono');
		for k=1:length(ix)
			if ix{k},
				f0(k) = [];
			end;	
		end;		
	end; 

	if 0,
	elseif length(f0)==1,
		feog = fullfile(p,f0.name);
		LAB_ID = 'BBCI'; 
	else
		f0 = dir(fullfile(p,'EOG/*.gdf')); % Graz BCI files 
		LAB_ID = 'GrazBCI'; 
		if length(f0)==1,
			feog = fullfile(p,'EOG',f0.name);
		else 	
			f0 = dir(fullfile(p,'*03.bkr')); % Graz BCI files 
			if length(f0)==0,
				f0 = dir(fullfile(p,'*12.bkr')); % Graz BCI files 
			end; 
			if length(f0)==1,
				feog = fullfile(p,f0.name);
				eegchan = 1:22;
			end
			LAB_ID = 'Graz22'; 
		end; 	
	end

	if length(f0)<1,
		feog = fn; 
		fprintf('error: no artefact file found. Use specified file.\n'); 
	end
	if length(f0)>1,
		fprintf('error: more than one artefact file found!\n'); 
		{f0.name}
		return; 
	end; 

	% load eye movement calibration data
        [s00,h0] = sload(feog);
	s0 = s00;

	FLAG.SEASON2DATA_PLAYER2 = 0;
	if strcmp(LAB_ID,'BBCI'),
		if ~isempty(strfind(Mode,'x'))	%%%% this is a hack for the SEASON2DATASET
			FLAG.SEASON2DATA_PLAYER2 = 1;
		end;

		% resting EEG with eyes open 
		ix = find([h0.EVENT.TYP==hex2dec('0114')] | [h0.EVENT.TYP==hex2dec('0115')]);
		r00 = [];
		for k=1:length(ix),
			r00 = [r00; s00(h0.EVENT.POS(ix(k))+(0:h0.EVENT.DUR(ix(k))),:)];
		end; 
			        
				
	        % find eye movements in BBCI recordings
        	tmp = bitand(2^15-1,h0.EVENT.TYP);
		tmp = find((tmp > hex2dec('0430')) & (tmp<=hex2dec('0439')));
		ix1 = min(tmp); 
		ix2 = max(tmp); 

	        % extract segment with large eye movements/EOG artifacts
	        s00 = s00(h0.EVENT.POS(ix1):h0.EVENT.POS(ix2)+h0.EVENT.DUR(ix2),:); 

	        if size(s00,1)<h0.SampleRate*10,
	        	fprintf(2,'WARNING GET_REGRESS_EOG: in %s no eye movements identified\n',h0.FileName);
	        end;
	        
	elseif strcmp(LAB_ID,'Graz22'),
		f0 = dir(fullfile(p,'*02.bkr')); % Graz BCI files 
		if length(f0)==0,
			f0 = dir(fullfile(p,'*11.bkr')); % Graz BCI files 
		end; 
		if length(f0)==1,
			feog = fullfile(p,f0.name);
		end
	        [r00,h0r] = sload(feog);

	else 
		r00 = [];         
	end;
	FLAG.HURST = ~isempty(strfind(Mode,'HURST')); % Vorobyov and Cichocki (2002)
	FLAG.AFOP = ~isempty(strfind(Mode,'AFOP')); % Boudet et al 2007
	FLAG.FOP  = ~FLAG.AFOP & ~isempty(strfind(Mode,'FOP')); % Boudet et al 2007
	FLAG.TDSEP1 = ~isempty(strfind(Mode,'TDSEP1'));
	FLAG.TDSEP2 = ~isempty(strfind(Mode,'TDSEP2'));
	FLAG.FFDIAG = ~isempty(strfind(Mode,'FFDIAG'));
	FLAG.Joyce = ~isempty(strfind(Mode,'JOYCE'));
	FLAG.Barbati = ~isempty(strfind(Mode,'BARBATI'));
	FLAG.Meinecke = ~isempty(strfind(Mode,'MEINECKE'));

	ix = strfind(Mode,'NGCA-');
	FLAG.NGCA = ~isempty(ix);
	if FLAG.NGCA,
		FLAG.NGCA = str2double(strtok(Mode(ix+5:end),' '));
	end;	

	ix = strfind(Mode,'TDSEP3-');
	FLAG.TDSEP3 = ~isempty(ix);
	if FLAG.TDSEP3,
		FLAG.TDSEP3 = str2double(strtok(Mode(ix+7:end),' '));
	end;	

	ix = strfind(Mode,'TDSEP-');
	FLAG.TDSEP = ~isempty(ix);
	if FLAG.TDSEP,
		FLAG.TDSEP = str2double(strtok(Mode(ix+6:end),' '));
	end;	

	ix = strfind(Mode,'ICA-');
	FLAG.ICA = ~isempty(ix);
	if FLAG.ICA,
		FLAG.ICA = str2double(strtok(Mode(ix+4:end),' '));
	end;	
	ix = strfind(Mode,'PCA-');
	FLAG.PCA = ~isempty(ix);
	if FLAG.PCA,
		FLAG.PCA = str2double(strtok(Mode(ix+4:end),' '));
	end;	
	if ~isempty(strfind(upper(Mode),'MSEC'));
		FLAG.PCA = 3;
	end;
	
	FLAG.REG = ~isempty(strfind(upper(Mode),'REG'));
        
	% downsampling 
	ix1 = strfind(upper(Mode),'FS=');
	if ~isempty(ix1),
		ix2 = strfind(upper(Mode(ix1:end)),'HZ');
		Fs  = Mode(ix1+3:ix1+ix2-2);
		Fs  = str2double(Fs);
		if length(Fs)==1 & isfinite(Fs),
			s00 = rs(s00,h0.SampleRate,Fs);
			h0.SampleRate = Fs; 
		else
			fprintf(2,'Warning GET_REGRESS_EOG: Fs="%s" could not be decoded. No filter is applied\n',tmp);
		end;	
	end;

	% filter data  
	ix1 = strfind(Mode,'FILT');
	if ~isempty(ix1),
		ix2 = strfind(upper(Mode(ix1:end)),'HZ');
		tmp = Mode(ix1+4:ix1+ix2-2);
		tmp((tmp=='-')|(tmp=='='))=' ';
		B   = str2double(tmp);
		if (length(B)==2) & all(isfinite(B))
			FLAG.Filter = B; 
			tmp = center(s00,1); tmp(isnan(tmp)) = 0;
		        tmp = fft(tmp); f=[0:size(s00,1)-1]*h0.SampleRate/size(s00,1);
		        w   = ((f>B(1)) & (f<B(2))) | ((f>h0.SampleRate-B(2)) & (f<h0.SampleRate-B(1)));
		        tmp(~w,:) = 0;
		        s00 = real(ifft(tmp));
		        clear tmp; 
		else
			fprintf(2,'Warning GET_REGRESS_EOG: Filter="%s" could not be decoded. No filter is applied\n',tmp);
		end;			
	end;
		
	%%%%% define EEG and EOG channels %%%%%
	CHANTYP = repmat(' ',h0.NS,1);
	CHANTYP([(h0.LeadIdCode>=996) & (h0.LeadIdCode<=1302)])='E'; % EEG
	if strcmp(LAB_ID, 'Graz22');
		CHANTYP(1:22)='E';
	end;	
 	eogchan = identify_eog_channels(h0);
        CHANTYP(any(eogchan,2)) = 'O';	% EOG
      	if FLAG.SEASON2DATA_PLAYER2,	
	        CHANTYP = CHANTYP([65:128,1:64]);
	        eogchan = eogchan([65:128,1:64],:);
	end;

	if isempty(strfind(upper(Mode),'MSEC1'))
		% MSEC,MSEC2 [default]
		CHAN = find((CHANTYP=='E') | (CHANTYP=='O'));
	else  % MSEC1 
		CHAN = find(CHANTYP=='E');
	end;
	chan = CHAN;
        nc   = length(CHAN); 

	%%%%% identify EOG components %%%%%	
	if 0,
	
	elseif FLAG.AFOP,
		error('Mode AFOP not supported');
	
	elseif FLAG.FOP,
		% Optimal Projection (Boudet et al, 2007)
		chan = find(CHANTYP=='E' | CHANTYP=='O');

	        RES.FOP.w = [];  	
	        [P,D] = eig(covm(r00(:,chan),'D'),covm(s00(:,chan),'D'));
	        TH = 0.3; % Boudet et al. 2007, p 1989. 
	        eogchan = sparse(chan,1:length(chan),1,h0.NS,length(chan))*P(:,diag(D)<TH);
	        if size(eogchan,2)>5,
	        	eogchan=eogchan(:,1:5); 
	        end;
	        RES.FOP.w = eogchan;  	

	elseif ~isempty(strfind(Mode,'bf-'))	
		% beamformer approach
		T = load('EOG_DIPOLS3'); 	% load lead field matrix
		T.Label = T.sa.clab_electrodes;
		T = leadidcodexyz(T);
		V = zeros(h0.NS,3);
		for k=find(h0.LeadIdCode)',
			ix = find(T.LeadIdCode==h0.LeadIdCode(k));
			if length(ix),
				V(k,:) = T.v(ix,4:6);
			end;	
		end;	
		eogchan = pinv(V)';
		if FLAG.SEASON2DATA_PLAYER2,
			eogchan = eogchan([h0.NS/2+1:h0.NS,1:h0.NS/2],:);
		end;
		
	elseif ~isempty(strfind(Mode,'bf+'))	
		% beamformer approach
		T = load('EOG_DIPOLS3');  % load lead field matrix
		T.Label = T.sa.clab_electrodes;
		T = leadidcodexyz(T);
		V = zeros(h0.NS,3);
		for k=find(h0.LeadIdCode)',
			ix = find(T.LeadIdCode==h0.LeadIdCode(k));
			if length(ix),
				V(k,:) = T.v(ix,4:6);
			end;	
		end;	
		R = eye(64);
		R(1:54,1:54)=R(1:54,1:54)-1/54;
		eogchan = [R*pinv(V(1:64,:))';zeros(64,3)];
		if FLAG.SEASON2DATA_PLAYER2,
			eogchan = eogchan([h0.NS/2+1:h0.NS,1:h0.NS/2],:);
		end;

	elseif FLAG.Meinecke,
		CHAN    = find(CHANTYP=='E');
	        chan    = sparse(CHAN,1:length(CHAN),1,h0.NS,length(CHAN));
		r       = [chan,eogchan];  % CAR and bipolar EOG
		tmp = s00*r; 
		tmp = tmp(~any(isnan(tmp),2),:); 

%     opt.verbose: 0 = no progress reports
%                  1 = progress reports (default)
%                  2 = progress reports and plot
%     opt.title:   string containing title of the data set
%     opt.B:       Bootstrap size (default B=100)
%     opt.sel:     time lag values for TDSEP

		opt.verbose = 0; 
		opt.title = 'BSS resampling';
		opt.B = 10; 
		opt.sel = 0:round(.15*h0.SampleRate)
		if ~isempty(strfind(upper(Mode),'JADE'));
			alg = 'jade';
		else 
			alg = 'tdsep';
		end;

		ica = ica_resamp(tmp',alg,opt)

%     ica.x:    the given mixtures
%     ica.A:    the estimated mixing matrix
%     ica.s:    the estimated source signals
%     ica.S:    the separability matrix
%     ica.perm: permutation of the estimated sources that makes the
%               block structure visible
		
		
		ochan = pinv(ica.A)';
		eogchan = r*ochan;


	elseif ~isempty(strfind(Mode,'JOYCE'))
		CHAN    = find(CHANTYP=='E');
	        chan    = sparse(CHAN,1:length(CHAN),1,h0.NS,length(CHAN));
		r       = [chan,eogchan];  % mono and bipolar EOG
		tmp = s00(~any(isnan(s00),2),:);

		sel = round(.15*h0.SampleRate);
		meth = 'TDSEP0';
		W1 = bss(zscore(tmp*[chan, eogchan]),meth,[],sel);
		W2 = bss(zscore(tmp*[chan,-eogchan]),meth,[],sel);

		% step 2
		c1 = corrcoef(tmp*[chan, eogchan]*W1,tmp*[chan, -eogchan]*W2);
		[sel2,ix]=find(c1 < -0.5);

		% step 3
		c2 = corrcoef(s00*eogchan,s00*[chan, eogchan]*W1);
		sel3 = find(any(abs(c2)>.3));
		
		% step 4
		c3 = rms(diff(zscore(s00*[chan, eogchan])*W1(:,sel3),[],1));
		sel4 = find(c3 < 0.2 * 500 / h0.SampleRate); 	
		%sel2',sel3,sel4,
		sel = unique([sel2(:)',sel3(sel4)]);
		FLAG.Joyce.C1 = c1;
		FLAG.Joyce.C2 = c2;
		FLAG.Joyce.C3 = c3;

		eogchan = [chan, eogchan]*W1(:,sel); 

	elseif ~isempty(strfind(Mode,'KIERKELS'))
		[t,r]=strtok(Mode);
		[meth]=strtok(r);
		
		eegchan = find(CHANTYP=='E');
		eogchan = find(CHANTYP=='O');
		chan    = [eegchan(:);eogchan(:)];
		tmp = s00(:,chan);
		tmp = tmp(~any(isnan(tmp),2),:);

		sel = round([0,1,2,3,5,10,20]/250*h0.SampleRate);
		W   = bss(tmp,meth,[],sel);
		A   = inv(W);
		r   = corrcoef(s00(:,chan)*W,s00(:,eogchan));
		sel = any(r>0.7, 2);
		eogchan = sparse(chan,1:length(chan),1,h0.NS,length(chan)) * W * A(:,sel);
		
		
	elseif ~isempty(strfind(Mode,'TING'))
		eegchan = find(CHANTYP=='E');
		eogchan = find(CHANTYP=='O');
		chan    = [eegchan(:);eogchan(:)];
		tmp = s00(:,chan);
		tmp = tmp(~any(isnan(tmp),2),:);

		sel = round(.15*h0.SampleRate);
		meth= 'TDSEP0';
		W  = bss(tmp,meth,[],sel);
		A  = inv(W);
		
		c1 = max(abs(tmp),[],1) * max(abs(A),[],2);  
		up = chan(length(eegchan)+1);
		lo = chan(length(eegchan)+2);
		le = chan(end-1);
		ri = chan(end);
		c2 = abs(A(:,up)-A(:,lo)) - max(max(A(:,1:length(eegchan))));
		c3 = abs(A(:,le)-A(:,ri)) - max(max(A(:,1:length(eegchan))));
		
		sel=zeros(1,3); 
		[tmp,sel(1)]=max(c1);
		[tmp,sel(2)]=max(c2);
		[tmp,sel(3)]=max(c3);
		sel=unique(sel);
		%find(c1>100),find(c2>0),find(c3>0),
		sel=[find(c1>100),find(c2>0),find(c3>0)]
		eogchan = sparse(chan,1:length(chan),1,h0.NS,length(chan)) * W * A(:,sel);
		

	elseif (FLAG.PCA>0) | (FLAG.ICA>0) | (FLAG.NGCA>0) | (FLAG.TDSEP>0) |(FLAG.TDSEP1>0) | (FLAG.TDSEP2>0) | (FLAG.TDSEP3>0) | (FLAG.HURST) | (FLAG.FFDIAG),
		% identify channels used for PCA

		if FLAG.PCA,
		        chan    = sparse(chan,1:nc,1,h0.NS,nc)*(eye(nc) - 1/nc); % Common Average Reference (CAR)
			r       = [chan,eogchan];  % CAR and bipolar EOG
			tmp = s00*r; 
			tmp = tmp(~any(isnan(tmp),2),:); 
			
			[u,s,ochan] = svds(tmp, FLAG.PCA); % get k (i.e. FLAG.PCA) largests PCA's
			eogchan = r*ochan;
			RES.PCA.w = eogchan; 
			
		elseif FLAG.ICA,
		        chan    = sparse(chan,1:nc,1,h0.NS,nc)*(eye(nc) - 1/nc); % Common Average Reference (CAR)
			r       = [chan,eogchan];  % CAR and bipolar EOG
			tmp = s00*r; 
			tmp = tmp(~any(isnan(tmp),2),:); 

			[A] = jade(tmp', FLAG.ICA); % get k (i.e. FLAG.PCA) largests PCA's
			ochan = pinv(A)';
			eogchan = r*ochan;
			
		elseif FLAG.HURST,
		        chan    = sparse(chan,1:nc,1,h0.NS,nc)*(eye(nc) - 1/nc); % Common Average Reference (CAR)
			r       = [chan,eogchan];  % CAR and bipolar EOG
			tmp = s00*r; 
			tmp = tmp(~any(isnan(tmp),2),:); 

			A  = jade(tmp', min(20,size(tmp,1))); 
			U  = pinv(A)';
			H  = hurst(s0*r*U)*log(2)  % account for c=1/2
			ix = find((H>=.58) & (H<=0.64));
			eogchan  = r*U(:,ix);
			h0.HURST = H;

		elseif FLAG.NGCA,	% NonGaussian Component Analysis
			% full rank required
			par.nbng = FLAG.NGCA;
		        chan = sparse(chan,1:nc,1,h0.NS,nc); 
			tmp = s00*chan; 
			tmp = tmp(~any(isnan(tmp),2),:); 
			[ngmatrix, projdata, signalmatrix] = NGCA(tmp',par);
			eogchan = chan*ngmatrix;

		elseif FLAG.FFDIAG,	% 
		        chan = sparse(chan,1:nc,1,h0.NS,nc); 
			tmp = s00*chan; 
			tmp = tmp(~any(isnan(tmp),2),:); 
			
			sel = 0:round(.15*h0.SampleRate);
			C= zeros(nc,nc,length(sel));
			for k= 1:length(sel)
			    C(:,:,k)= tdCorr(tmp', sel(k));
			end
			[W,d] = ffdiag2(C);
			W = (W./repmat(sqrt(sum(W.*W,2)),1,nc))';
			[tmp,ix]=sort(-rms(s00*chan*W,1));
			eogchan = chan*W(:,ix(1:3));

		elseif FLAG.TDSEP,	% 
		        chan = sparse(chan,1:nc,1,h0.NS,nc)*(eye(nc) - 1/nc); % Common Average Reference (CAR)
			tmp = s00*chan; 
			tmp = tmp(~any(isnan(tmp),2),:); 

			sel = round([0,1,2,3,5,10,20]/256*h0.SampleRate);
			%%%  DO NOT USE
			[W,d] = tdsep(tmp',sel);
			eogchan = chan*W(1:3,:)';

		elseif FLAG.BSS,	% 
		        chan1= sparse(chan,1:nc,1,h0.NS,nc);
			tmp = s00*chan1; 
			tmp = tmp(~any(isnan(tmp),2),:); 
			sel = 0:round(.15*h0.SampleRate);
			W   = bss(tmp',alg,[],sel);

			%%% select 1st 3 components  %%%
			eogchan = chan*W(1:3,:)';

			%%% select 3 components with minumum angle to b0 %%%
			[h0.REGRESS, s01] = regress_eog(s00,chan,eogchan); 
			iW  = inv(W); 
			for k = 1:size(W,2),
				phi(k) = subspace(h0.REGRESS.b0(2:3,:)',iW(:,k));
			end; 	
			[phi,ix] = sort(phi); 
			eogchan = chan1*W(ix(1:3),:)';
			
		elseif FLAG.TDSEP1,	% 
		        chan1= sparse(chan,1:nc,1,h0.NS,nc);
			tmp = s00*chan1; 
			tmp = tmp(~any(isnan(tmp),2),:); 
			sel = 0:round(.15*h0.SampleRate);
			W   = tdsep0(tmp',sel);
			iW  = inv(W); 
			[h0.REGRESS, s01] = regress_eog(s00,chan,eogchan); 

			for k = 1:size(W,2),
				phi(k) = subspace(h0.REGRESS.b0(2:3,:)',iW(:,k));
			end; 	
			[phi,ix] = sort(phi); 
			eogchan = chan1*W(ix(1:3),:)';
			
		elseif FLAG.TDSEP2,	% 
		        chan1= sparse(chan,1:nc,1,h0.NS,nc);
			tmp = s00*chan1; 
			tmp = tmp(~any(isnan(tmp),2),:); 
			sel = 0:round(.15*h0.SampleRate);
			W   = tdsep0(tmp',sel);
			iW  = inv(W); 
			[h0.REGRESS, s01] = regress_eog(s00,chan,eogchan); 

			c  = [];
			ix = [];
			P = eye(length(chan)); 
			for k0 = 1:3, 
				for k = 1:size(W,2),
					phi(k) = subspace(P*h0.REGRESS.b0(2:3,:)',P*iW(:,k));
				end; 	
				[sp,ix(k0)] = min(phi);
				w = iW(:,ix);
				c = [c,iW(:,ix)]; 
				P = P*[eye(nc)-w*inv(w'*w)*w'];
			end; 	 
			eogchan = chan1*W(ix,:)';
			
		elseif FLAG.TDSEP3,	% 
			% full rank required
		        chan    = sparse(chan,1:nc,1,h0.NS,nc); 
			r       = [chan,eogchan];  % CAR and bipolar EOG
			tmp = s00*chan; 
			tmp = tmp(~any(isnan(tmp),2),:); 

			sel = round([0,1,2,3,5,10,20]/256*h0.SampleRate);
			[y,W,d] = tdsep3(tmp',sel);
			eogchan = chan*W(1:3,:)';
		end;
		 
		 
	elseif FLAG.REG, %REGRESSION
	        % define xEOG channels
	        FLAG.REG = 1; 
		if ~isempty(strfind(upper(Mode),'REG+CAR1'));
	        	eogchan = [eogchan,any(eogchan,2)];
		elseif ~isempty(strfind(upper(Mode),'REG+CAR2'));
			eogchan = [eogchan,double([CHANTYP=='O']/sum(CHANTYP=='O')-[CHANTYP>' ']/sum(CHANTYP>' '))];
		end;
	end;        

       	% remove EOG channels for list of corrected channels
	% chan = 1:h0.NS;
	chan = find(CHANTYP=='E');
        % regression analysis - compute correction coefficients. 
	[h0.REGRESS, s01] = regress_eog(s00,chan,eogchan); 
	h0.REGRESS.FLAG = FLAG;	

	%%%%% post-regression improvement
	if length(Mode)>6,
	if ~isempty(strfind(Mode([1:4,6,7]),'REG+CA'))	
		% select EEG and EOG channels and do CAR
		echan = find(CHANTYP>' ');
	        nc    = length(echan); 
	        echan = sparse(echan,1:nc,1,h0.NS,nc)*(eye(nc) - 1/nc); % Common Average Reference (CAR)

		tmp = strtok(Mode);
		nc = str2double(tmp(8:end));
		if isempty(nc) | isnan(nc),
			nc = 1;
		else
		end;		
        	tmp = center(s01,1)*echan;
		tmp = tmp(~any(isnan(tmp),2),:); 
		if ~isempty(strfind(Mode,'REG+PCA'))	
			% get largest PCA's
			[u,s,zeog] = svds(tmp, nc);
	        elseif ~isempty(strfind(Mode,'REG+ICA'))
	        	%%% ICA: identify one component
			[A,s] = jade(tmp',nc);
			zeog  = pinv(A)';
		else 	
			fprintf(2,'Mode %s not supported - regression is used only.\n',Mode); 
			return;
		end; 

		w = sparse([eogchan,h0.REGRESS.r0*echan*zeog]);
%rank(full(w))
w(isnan(w))=0;
h0.REGRESS.r0,
%rank(full(w))
		[h0.REGRESS,s01] = regress_eog(s00,chan,w);
	end
	end
	h0.REGRESS.FLAG = FLAG;	


	%%% quantifying the noise 
