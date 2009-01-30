function H=plota(X,arg2,arg3,arg4,arg5,arg6,arg7)
% PLOTA plots all kind of data types
%
% PLOTA(Labels,'ELPOS')
%       displays 2-D location of EEG channels described by Labels
% PLOTA(Labels,'ELPOS3')
%       displays 3-D location of EEG channels described by Labels
%
% PLOTA(X [,Mode])
%
% if X.TYPE=='EVENT' and X.EVENT
% 
%
% X.datatype determines type of data
%    DATATYPE   Mode
%   'MVAR'      'AutoSpectrum'
%   'MVAR'      'SPECTRUM','logS'
%   'MVAR'      'logH'
%   'MVAR'      'Phase'
%   'MVAR',	'COHERENCE'
%   'MVAR'      'iCOH'	imaginary coherence
%   'MVAR'      'iSpectrum'	imaginary part of spectrum
%   'MVAR'      'DTF'
%   'MVAR'      'PDC'
%   'MVAR'      'dDTF'
%   'MVAR'      'ffDTF'
%   'TF-MVAR'	Time-frequency MVAR analysis
%		e.g. plota(X, 'PDC', alpha [, Y]);	%
%
%   'MEAN+STD'
%       plota(X,hf,minmean,maxmean,maxstd [,trigger])
%       arg1 ... R	
%       arg2 ... hf (handles to figures)
%       arg3 ... minmean (minimum of mean)
%       arg4 ... maxmean (maximum of mean)
%       arg5 ... maxstd (maximum of standard deviation)
%       arg6 ... trigger (trigger instant) [optional]
%
%   'QRS_events' shows inter-beat-interval (IBI) of the ecg 
%
%   'HISTOGRAM'	'log'	chansel
%   'HISTOGRAM'	'log+'	chansel
%   'HISTOGRAM'	'log '	chansel
%   'HISTOGRAM'	'lin'	chansel
%
%   'SIESTA_HISTOGRAM'	chansel
%
%   'DBI-EPs'
%   'TSD1'
%   'TSD_BCI7'
%   'MDist-matrix'
%   'MD'
%   'SCATTER'
%   'STAT2'
%   ''
%   ''
%   'REV' Mode='3D'
%   'REV' Mode='2D'
%
% REFERENCE(S):


%	$Id: plota.m,v 1.1 2009-01-30 06:04:51 arno Exp $
%	Copyright (C) 2006,2007,2008 by Alois Schloegl <a.schloegl@ieee.org>
%       This is part of the BIOSIG-toolbox http://biosig.sf.net/
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


h = [];
if nargin>1,
        if strncmpi(arg2,'ELPOS',5),
                if ~isstruct(X),
                        tmp = X; X=[]; 
                        X.Label = cellstr(tmp);
                end;
                X.datatype = upper(arg2);
                clf;
        elseif strcmpi(arg2,'SCATTER'),
        	tmp = X; X=[];
        	X.datatype = upper(arg2);
        	if isnumeric(tmp)
	        	X.data  = tmp;
	        	if nargin<3, arg3='x'; end;
	        	plota(X,arg3);
	        end;
	        	
        end;
end;

if isfield(X,'datatype');

elseif isfield(X,'TYPE');
        if strcmp(X.TYPE,'EVENT')
                if any(X.EVENT.TYP==hex2dec('0501'));
                        X.datatype = 'QRS_events';
                        H = X;
                end
                
        elseif strcmp(X.TYPE,'ELPOS')
                if isfield(X,'ELEC') & isfield(X,'Label');
                        if isfield(X.ELEC,'XYZ');
                                X.datatype = 'ELPOS'; 
                        end;
                end;
        end;
else
        return;
end;

if 0,

elseif strcmp(X.datatype,'MVAR'),
        if ~isfield(X,'A') || ~isfield(X,'B'),
                fprintf(2,'Error PLOTA: MVAR missing input data\n');
                return;
        end;

        [K1,K2] = size(X.A);
        p = K2/K1-1;
        [K1,K2] = size(X.B);
        q = K2/K1-1;
        if ~isfield(X,'C');
                X.C=ones(K1,K1);
        end;

        if isfield(X,'SampleRate');
                Fs = X.SampleRate;
        elseif nargin < 4,
                Fs = 1; pi*2;
        else
                Fs = arg4;
        end;

        if nargin<2,
                Mode= 'DTF';
                Fs  = 1;
        else
                Mode = arg2;
        end;

        if nargin<3,
                N=512;
        else
                N=arg3;
        end;

        if all(size(N)==1)
                f = (0:N)/(2*N)*Fs;
        else
                f = N;
                N = length(N);
        end;

        if isfield(X,'Label');
                Label = cellstr(X.Label);
        elseif isfield(X,'ElectrodeName');
                if isstruct(X.ElectrodeName);
                        Label = X.ElectrodeName;
                else
                        for k=1:K1,
                                Label{k}=X.ElectrodeName(k,:);
                        end;
                end;
        else
                for k=1:K1,
                        Label{k}=sprintf('#%02i',k);
                end;
        end;

        if strcmpi(Mode,'Eigen'),
		[S, Serr, per, tau, exctn, lambda] = armode2(-X.A(:,K1+1:end), X.C);
		[per(1,:)',Fs./per(1,:)',tau(1,:)'/Fs,exctn',lambda]
		return;
	end; 

        [S,h,PDC,COH,DTF,DC,pCOH,dDTF,ffDTF, pCOH2, PDCF, coh, GGC, Af, GPDC, GGC2]=mvfreqz(X.B,X.A,X.C,f,Fs);
        Phase = zeros(size(h));
        for k1=1:K1;
        for k2=1:K2;
        	Phase(k1,k2,:) = unwrap(squeeze(angle(S(k1,k2,:))))*180/pi;
        end;
        end;
	dT = Phase./reshape(repmat(360*f(:)',[K1*K1,1]),[K1,K1,length(f)]);

        range = [0,1]; % default range
        if ~isempty(strfind(Mode,'Auto')),
                if 0,

                elseif strcmpi(Mode,'AutoSpectrum'),
                        R = abs(S);
                        range = [min(R(:)),max(R(:))];
                        range = [1e-8,1e0];
                else
                        error(['unknown MVAR-parameter: ',arg2])
                end;
                M = size(S,1);
                K1 = ceil(sqrt(M));
                K2 = ceil(M/K1);
                for k1=1:M;
                        subplot(K1,K2,k1);
                        semilogy(f,squeeze(R(k1,k1,:)));
                        axis([0,max(f),range]);
                        ylabel(Label{k1});
                end;

        else
                if 0,

                elseif strcmpi(Mode,'Spectrum') || strcmpi(Mode,'logS'),
                        R = abs(S);
                        range = [min(R(:)),max(R(:))];
                        range(1) = min(range(1),range(2)/100);
                        %range=[1e-2,1e2];
                elseif strcmpi(Mode,'logh'),
                        R = abs(h);
                        range = [min(R(:)),max(R(:))];
                        range(1) = min(range(1),range(2)/100);
                        %range=[1e-2,1e2];
                elseif strcmpi(Mode,'iSpectrum'),
                        R = imag(S);
                        range = [min(R(:)),max(R(:))];
                elseif strcmpi(Mode,'rSpectrum'),
                        R = real(S);
                        range = [min(R(:)),max(R(:))];
                elseif strcmpi(Mode,'Phase') | strcmpi(Mode,'phaseS') ,
                        R = zeros(size(S));
                        for k1=1:K1;
                                for k2=1:K2;
                                        R(k1,k2,:) = unwrap(squeeze(angle(S(k1,k2,:))))*180/pi;
                                end;
                        end;
                        range = [-180,180]*2;
                        range = [min(R(:)),max(R(:))];
                elseif strcmpi(Mode,'Phase') | strcmpi(Mode,'phaseh') ,
                        R = zeros(size(h));
                        for k1=1:K1;
                                for k2=1:K2;
                                        R(k1,k2,:) = unwrap(squeeze(angle(h(k1,k2,:))))*180/pi;
                                end;
                        end;
                        range = [-180,180]*2;
                        range = [min(R(:)),max(R(:))];
                elseif strcmpi(Mode,'PDC'),
                        R = PDC;
                elseif strcmpi(Mode,'GPDC'),
                        R = GPDC;
                elseif strcmpi(Mode,'Coherence') | strcmpi(Mode,'COH'),
                        R = abs(COH);
                elseif strcmpi(Mode,'iCOH') | strcmpi(Mode,'imagCOH'),
                        R = imag(COH);
                        range = [-1,1];
                elseif strcmpi(Mode,'pCOH'),
                        R = abs(pCOH);
                elseif strcmpi(Mode,'pCOH2'),
                        R = abs(pCOH2);
                elseif strcmpi(Mode,'coh'),
                        R = abs(coh);
                elseif strcmpi(Mode,'icoh'),
                        R = imag(coh);
                        range = [-1,1];
                elseif strcmpi(Mode,'GGC'),
                        R = log10(GGC);
                        range = [min(R(:)),max(R(:))];
                        range = [.1,max(R(:))];
                elseif strcmpi(Mode,'GGC2'),
                        R = (GGC2);
                        range = [min(R(:)),max(R(:))]
                        %range = [.1,max(R(:))];
                elseif strcmpi(Mode,'Af'),
                        R = abs(Af);
                        for k=1:size(R,1),
%                        	R(k,k,:)=NaN;
                        end; 	
                        range = [min(R(:)),max(R(:))].*[.9,2]
                        %range = [[.001,1]*max(R(:))]
                        range = [0,max(R(:))];
                        R = abs(Af);
                elseif strcmpi(Mode,'Af1'),
                        R = log10(abs(Af))+3;
                        range = [min(R(:)),max(R(:))]
                        %range = [[.001,1]*max(R(:))]
                elseif strcmpi(Mode,'PDCF'),
                        R = PDCF;
                elseif strcmpi(Mode,'DTF'),
                        R = DTF;
                elseif strcmpi(Mode,'dDTF'),
                        R = dDTF;
                elseif strcmpi(Mode,'ffDTF'),
                        R = ffDTF;
                elseif strcmpi(Mode,'dT'),
                        R = dT;
                        tmp = dT(isfinite(dT(:)));
                        range = [min(tmp(:)),max(tmp(:))]
                elseif strcmpi(Mode,'DCF1'),
                        R = S;
                        for k1=1:K1,
                                for k2=1:K1,
                                        R(k1,k2,:) = sqrt(X.C(k2,k2)/(2*pi*Fs))*abs(h(k1,k2,:))./sqrt(S(k1,k1,:));
                                end;
                        end;
                        range = [0,1];
                        %			range = [min(R(:)),max(R(:))];
                elseif strcmpi(Mode,'DCF2'),
                        R = S;
                        for k1=1:K1,
                                for k2=1:K1,
                                        R(k1,k2,:) = abs(h(k1,k2,:))./sqrt(S(k1,k1,:));
                                end;
                        end;
                        range = [0,1];
                elseif strcmpi(Mode,'DCF'),
                        tmp = S;
                        for k=1:K1,
                                tmp(k,:,:) = tmp(k,repmat(k,1,K1),:);
                        end;
                        R = h./sqrt(tmp);
                        range = [-1,1];
                        %			range = [min(R(:)),max(R(:))];
                else
                        error(['unknown MVAR-parameter: ',arg2])
                end;

                for k1=1:K1;
                        for k2=1:K2;
                                subplot(K1,K2,k2+(k1-1)*K1);
                                if strcmpi(Mode,'logS') %| strcmpi(Mode,'Af'),
                                        semilogy(f,squeeze(R(k1,k2,:)));
                                else
                                        area(f,squeeze(R(k1,k2,:)));
                                end;
                                axis([0,max(f),range]);
                                if k2==1;
                                        ylabel(Label{k1});
                                end;
                                if k1==1;
                                        title(Label{k2});
                                end;
                        end;
                end;
        end;
        if exist('suptitle','file')
%                suptitle(Mode);
        end;


elseif strncmp(X.datatype,'TF-MVAR',7)    % logS2 and S1

        %GF = {'C','DC','AR','PDC','DTF','dDTF','ffDTF','COH','pCOH','pCOH2','S','h','phaseS','phaseh','coh','logh','logS'};

        if nargin<2,
                arg2 = 'logS1';
        end;
        if nargin<3,
                alpha = .01;
        elseif isnumeric(arg3),
                alpha = arg3;
        elseif isempty(str2num(arg3))
                alpha = flag_implicit_significance;
        else
                alpha = arg3;
        end;
        if nargin<4,
                Y = [];
        else
                Y = arg4;
        end;
        if isfield(X,'SampleRate') & (size(X.T,1)>1),
		XT = (X.T(1,:)+X.T(2,:))/(2*X.SampleRate);
	else
		XT = X.T; 
	end;
        if isfield(Y,'SampleRate') 
        if (size(Y.T,1) > 1),
		YT = (Y.T(1,:)+Y.T(2,:))/(2*Y.SampleRate);
	end;
	end;
	
        gf = arg2;
        AUTO = ~isempty(strfind(lower(gf),'auto'));
        if AUTO,
                gf = strrep(gf,'Auto','');
        end;
        GF = strtok(gf);
	if strcmpi(arg2,'Eigen'),
		[K1,K2] = size(X.M.AR);
		Fs = X.SampleRate;
		[S, Serr, per, tau, exctn, lambda] = armode2(-X.M.AR(:,K1+1:end), X.M.C);
		fprintf(1,'Periode [s]\tf0 [Hz]    \ttau [s] \texcitation [?]\tlambda\n')
		fprintf(1,'%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\n',[per(1,:)',Fs./per(1,:)',tau(1,:)'/Fs,exctn',lambda]');
		return;	
	elseif ~isfield(X.M,GF)
                warning('PLOTA TFMVAR_ALL: field %s is unknown\n',GF);
        end;
        MONO = strcmp(GF,'logS1') | strcmp(GF,'S1');
        if strcmp(GF,'AR1') | strcmp(GF,'C1');
                MONO=1;
                f = 1:size(getfield(X.M,GF),2);
        else
                f = X.F;
        end;

        fidx = repmat(logical(1),size(f)); 	% selection of frequency band
        nix  = repmat(logical(1),size(XT));
        if 0,
                % ERD MVAR chapter
                warning('use project-specific setting: f<45Hz, t<7.0s')
                fidx = (f<=45)&(f>0); 	% selection of frequency band
                nix(XT>=7) = 0; 	% selection of segments (in time)
                nix(XT< 2.0) = 0; 	% selection of segments (in time)
                nix(1) = 0; 		% hack to remove reference segment
                nix = logical(nix); 
        end;
        nix(1) = logical(0); 		% hack to remove reference segment
        nix = logical(nix); 

        M   = size(X.M.AR,1);
        tmp = size(X.M.AR);
        MOP = tmp(2)/tmp(1);

        if ~isfield(X,'Label'),
                for k1 = 1:M,
                        Label{k1}=['# ',int2str(k1)];
                end;
        elseif ischar(X.Label),
                Label = cellstr(X.Label);
        else
                Label = X.Label;
        end;
        nr = ceil(sqrt(M));
        nc = ceil(M/nr);
        if 0, nargin>2,
                hf = arg3;

        elseif AUTO | MONO,
                NR = ceil(sqrt(M));
                NC = ceil(M/NR);
                for k1 = 1:M,
                        hf(k1) = subplot(NR,NC,k1);
                end;

        elseif M*M<200,
                for k1 = 1:M,
                        for k2 = 1:M,
                                %				hf(k1,k2)=subplot(M,M,(k2-1)*M+k1);	%transposed
                        end;
                end;
        else
                subplot(111);
        end;

        if strcmp(X.datatype,'TF-MVAR')
                %		warning('TF-MVAR 1.0 has');
                Xnormfactor = X.N;
        else
                Xnormfactor = 1;
                Ynormfactor = 1;
        end;

	type = ''; 
        if ~isempty(strfind(gf,'eventrelated'))

                [tmp]  = str2double(gf);
                [gf,r] = strtok(gf);
                [t2,r] = strtok(r);
                [t3,r] = strtok(r);

        %suche das erste vorkommen des Zeitfenstermittelpunktes, der als
        %Parameter hinter 'eventrelated' übergeben wurde - dieses TW ist
        %das Referenzfenster für den eventrelated-Vergleich (meistens BL)
		ix  = min( find(XT == max(tmp) ));
%                nix(ix) = 0;
                rix = repmat(ix,1,sum(nix));

                m   = real(getfield(X.M,gf));
                se  = real(getfield(X.SE,gf))*Xnormfactor;
                if MONO,
                        x0  = m(:,fidx,nix) - m(:,fidx,rix);
                        ci0 = sqrt(se(:,fidx,nix).^2 + se(:,fidx,rix).^2);
                else
                        x0  = m(:,:,fidx,nix) - m(:,:,fidx,rix);
                        ci0 = sqrt(se(:,:,fidx,nix).^2 + se(:,:,fidx,rix).^2);
                end;
                x0(~isfinite(x0))=NaN;
                clim = [-1,1]*max(abs(x0(:)));
                if isfield(X,'SampleRate')
	                type = [' - eventrelated with : ' , sprintf('%i - %i s', X.T(:,ix)/X.SampleRate) ];
		else
	                type = [' - eventrelated with : ' , sprintf('%i s', X.T(ix)) ];
		end; 
		
        elseif isempty(Y);
                x0  = real(getfield(X.M,gf));
                ci0 = getfield(X.SE,gf)*Xnormfactor;
                if MONO,
                        x0  = x0(:,fidx,nix);
                        ci0 = ci0(:,fidx,nix);
                else
                        x0  = x0(:,:,fidx,nix);
                        ci0 = ci0(:,:,fidx,nix);
                end;
                x0(~isfinite(x0))=NaN;
                clim = [min(x0(:)),max(x0(:))];

        else
                x0  = (real(getfield(X.M,gf)) - real(getfield(Y.M,gf)));
                ci0 = sqrt((real(getfield(X.SE,gf))*Xnormfactor).^2 + (real(getfield(Y.SE,gf))*Ynormfactor).^2);

                x0  = x0(:,:,fidx,nix);
                ci0 = ci0(:,:,fidx,nix);

                x0(~isfinite(x0))=NaN;
                clim = [-1,1]*max(abs(x0(:)));
        end;
        XT = XT(:,nix);

        tmp = diff(XT);
        if isempty(tmp)
        	sz = size(X.M.A);
                X.A = [eye(sz(1)), -reshape(X.M.A,[sz(1),sz(2)*sz(3)])];
                X.B = eye(size(X.M.A,1));
                X.C = X.M.C;
                % X.SampleRate = 2*max(X.F); 
                X.datatype='MVAR'; 

                m = getfield(X.M,arg2);
                se = getfield(X.SE,arg2);
		if strmatch(upper(arg2),{'PDC','DTF','COH','ICOH','PCOH'})
			range = [0,.4];
		elseif strcmpi(arg2,'Af')
			tmp=m;
			for k=1:size(m,1),
				tmp(k,k,:)=NaN;
			end; 	
	                range=[min(tmp(:)-se(:)),max(tmp(:)+se(:))]
		else		
	                range=[min(m(:)),max(m(:))];
	        end;         
                for k1=1:sz(1),
                for k2=1:sz(2),
                	subplot(sz(1),sz(2),(k1-1)*sz(2)+k2);
                	%area(X.F,squeeze(m(k1,k2,:)));
                	errorbar(X.F,squeeze(m(k1,k2,:)),squeeze(se(k1,k2,:)));
                	axis([min(X.F),max(X.F),range])
                	%set(gca,'YLIM',range);
                	if isfield(X,'Label')
	                	if k2==1, ylabel(X.Label{k1}); end; 
        	        	if k1==1, title(X.Label{k2}); end;
        	        end;
                end;
                end; 
                H = X; 
                fprintf(2,'Warning PLOTA: single time segment not implemented, yet.\n');
                return;
        elseif any(tmp-tmp(1))
                warning('time scale is not equidistant - xlabels are incorrect')
        end;

        if (strcmp(gf,'DC') | strcmp(gf,'C'))
                caxis(clim);
                cm = colormap;
                for k1 = 1:M,
                        for k2 = 1:M,
                                subplot(hf(k1*M-M+k2));
                                x  = x0(k1,k2,1,:);
                                ci = ci0(k1,k2,1,:);
                                if 0,
                                elseif alpha < .5,
                                        xc = 2 + round(62*(squeeze(x)-clim(1))/diff(clim));
                                        sz = size(x);
                                        %x = x(:);
                                        bf = prod(size(x));
                                        %xc(abs(x) < (ci*norminv(1-alpha/(2*bf)))) = 1;
                                        x(abs(x) < (ci*tinv(1-alpha/2,X.N))) = NaN;
                                        %x(abs(x) < .5) = NaN;
                                        %x = reshape(x,sz);
                                        cm(1,:) = [1,1,1];
                                        colormap(cm);
                                else
                                        xc = 1+round(63*(squeeze(x)-clim(1))/diff(clim));
                                        colormap('default');
                                end;

                                %h = semilogy(XT,[x(:),ci(:)]*[1,1,1;0,-1,1]);
                                h = plot(XT,[x(:),ci(:)]*[1,1,1;0,-1,1]);
                                set(h(1),'color',[0,0,1]);
                                set(h(2),'color',[0.5,0.5,1]);
                                set(h(3),'color',[0.5,0.5,1]);
                                v  = axis; v(1)=min(XT);v(2)=max(XT);axis(v);
                                axis([min(X.T),max(XT),clim])
                                %h  = imagesc(XT,X.F,squeeze(x),clim);
                                if k2==1, title(Label{k1}); end;
                                if k1==1, ylabel(Label{k2});end;
                        end;
                end;

        elseif AUTO || MONO,
                caxis(clim);
                cm = colormap;

                for k1 = 1:M,
                        subplot(hf(k1));
                        if MONO,
                                x  = x0(k1,:,:);
                                ci = ci0(k1,:,:);
                        else
                                x  = x0(k1,k1,1:length(X.F),:);
                                ci = ci0(k1,k1,1:length(X.F),:);
                        end;
                        if alpha < .5,
                                xc = 2 + round(62*(squeeze(x)-clim(1))/diff(clim));
                                sz = size(x);
                                %x = x(:);
                                bf = prod(size(x));
                                %xc(abs(x) < (ci*norminv(1-alpha/(2*bf)))) = 1;
                                xc(abs(x) < (ci*tinv(1-alpha/2,X.N))) = 1;
                                %x(abs(x) < .5) = NaN;
                                %x = reshape(x,sz);
                                cm(1,:) = [1,1,1];
                                colormap(cm);
                        else
                                xc = 1+round(63*(squeeze(x)-clim(1))/diff(clim));
                                colormap('default');
                        end;
                        x1 = reshape(cm(xc,1),size(xc));
                        x2 = reshape(cm(xc,2),size(xc));
                        x3 = reshape(cm(xc,3),size(xc));

                        %h = imagesc(XT,X.F,cat(3,x1,x2,x3)*diff(clim)+clim(1),clim);
                        h = imagesc(XT,f(fidx),cat(3,x1,x2,x3),clim);
                        %h  = imagesc(XT,X.F,squeeze(x),clim);
                        ylabel(Label{k1});
                end;

        elseif M*M<200,
                caxis(clim);
                cm = colormap;
                for k1 = 1:M,
                        for k2 = 1:M,
                                %				subplot(hf(k1*M-M+k2));		% display transposed because of predefined hf above
                                subplot(M,M,(k1-1)*M+k2);	% display transposed
                                %				subplot(M,M,(k2-1)*M+k1);
                                x  = x0(k1,k2,:,:);
                                ci = ci0(k1,k2,:,:);
                                if alpha < .5,
                                        xc = 2 + round(62*(squeeze(x)-clim(1))/diff(clim));
                                        sz = size(x);
                                        %x = x(:);
                                        bf = prod(size(x));
                                        xc(abs(x) < (ci*norminv(1-alpha/(2*bf)))) = 1;
                                        %xc(abs(x) <= (ci*tinv(1-alpha/2,X.N))) = 1;
                                        %x(abs(x) < .5) = NaN;
                                        %x = reshape(x,sz);
                                        cm(1,:) = [1,1,1];
                                        colormap(cm);
                                else
                                        xc = 1+round(63*(squeeze(x)-clim(1))/diff(clim));
                                        colormap('default');
                                end;
				xc(xc~=xc)=1;

                                x1 = reshape(cm(xc,1),size(xc));
                                x2 = reshape(cm(xc,2),size(xc));
                                x3 = reshape(cm(xc,3),size(xc));

                                %h = imagesc(XT,X.F,cat(3,x1,x2,x3)*diff(clim)+clim(1),clim);
                                h = imagesc(XT,X.F(fidx),cat(3,x1,x2,x3),clim);

                                if k1==1, title(Label{k2}); end;
                                if k1<M,
                                        set(gca,'xticklabel','');
                                end;
                                if k2==1,
                                        ylabel(Label{k1});
                                else
                                        set(gca,'yticklabel','');
                                end;
                                set(gca,'position',get(gca,'position').*[1,1,1.3,1.3]-[.05,0,0,0]);
                        end;
                end;
        else
                sz = size(x0);

                x  = reshape(permute(x0 ,[3,1,4,2]),[prod(sz([1,3])),prod(sz([2,4]))]);
                ci = reshape(permute(ci0,[3,1,4,2]),[prod(sz([1,3])),prod(sz([2,4]))]);
                % transposed
                %		x  = reshape(permute(x0 ,[3,2,4,1]),[prod(sz([2,3])),prod(sz([1,4]))]);
                %		ci = reshape(permute(ci0,[3,2,4,1]),[prod(sz([2,3])),prod(sz([1,4]))]);

		if ~exist('OCTAVE_VERSION','builtin');
	                caxis(clim);
	                cm = colormap;
	        else        
	                cm = colormap('default');
	        end;        


                if alpha < .5,
                	xc = 2 + round(62*(squeeze(x)-clim(1))/diff(clim));

                        %x = x(:);
                        bf = prod(size(x));
                        %xc(abs(x) < (ci*norminv(1-alpha/(2*bf)))) = 1;
                        xc(abs(x) <= (ci*tinv(1-alpha/2,X.N))) = 1;
                        %x(abs(x) < .5) = NaN;
                        cm(1,:) = [1,1,1];
                        colormap(cm);
                else
                        xc = 1+round(63*(squeeze(x)-clim(1))/diff(clim));
                        colormap('default');
                end;
		xc(xc~=xc)=1;
                x1 = reshape(cm(xc,1),size(xc));
                x2 = reshape(cm(xc,2),size(xc));
                x3 = reshape(cm(xc,3),size(xc));

                                
                % cat(3,x1,x2,x3) enthält das Bild, wobei es die Dimension
                % (m-by-n-by-3 besitzt --> Zeile * Spalte * RGB-Floats
                H.xKord = repmat(XT,1,sz(1)); %X-Koordinaten
                H.yKord = repmat(X.F(fidx),1,sz(2)); %Y-Koordinaten
                H.farbwerte = squeeze(cat(3,x1,x2,x3)); %Farbwerte                               

		if isfield(X,'Label')
	                H.Label = X.Label;
	        else  
	        	H.Label = cellstr(sprintf('#03i',[1:size(X.M.S,1)]'));
	        end; 	        
                H.type = type;
                H.method = gf;
                H.clim = clim;
                H.F = X.F;
                H.T = XT;
                H.zuPlottendeTW = nix;
                
                if(isfield(H,'cond'))                
                    H.cond = X.cond;
                else
                    H.cond = '';
                end;
                if(isfield(H,'subj'))
                    H.subj = X.subj;
                else
                    H.subj = '';
                end;

                                
% Method to plot the matrixes that were calculated with calc_TF_Matrix
% the input argument is the output of calc_TF_Matrix

if ~exist('OCTAVE_VERSION','builtin');
	imagesc(H.xKord, H.yKord, H.farbwerte, H.clim);               		
else
        cm = load('cm.txt');
        colormap(cm);
	imagesc(H.xKord*10, H.yKord, x);               		
end;	

%set data proporties--------------------------
numPlots = size(H.Label,1);   % Number of channels
numXIntervalle = size(H.T,2); % Anzahl der Intervalle auf der X-Achse je Kanal (Auflösung Zeitfenster)
numYIntervalle = size(H.F,2); % Anzahl der Intervalle auf der Y-Achse je Kanal (Frequenzen)


%Vertical lines----------------------------
axisTmp = axis; %Achsenminima und -maxima
xDim = axisTmp(2) - axisTmp(1);
xInterval = xDim / (numPlots); %Größe eine Subplots auf der X-Achse
xVerLines = [1:numPlots+1]  * xInterval + axisTmp(1); %X-Koordinaten der vertikalen Linien

v = get(gca);
for k=1:numPlots+1
    lh = line([xVerLines(k) xVerLines(k)] , [0 200]);
    set(lh,'Color',[.25 .25 .25]);
end;

%Horizontal lines-------------------------
yDim = axisTmp(4) - axisTmp(3);
yInterval = yDim / numPlots; %Größe eines Subplots auf der Y-Achse
HorLines = [1:numPlots]  * yInterval +  axisTmp(3); %Y-Koordinaten der horizontalen Linien

v = get(gca);
for k=1:numPlots
    lh = line([0 200],[HorLines(k) HorLines(k)]);
    set(lh,'Color',[.25 .25 .25]);
end;

%Y-Labels------------------------------------
  ylabel('') % Label ausschalten

%X-Labels------------------------------------
xlabel('')
title('')
for k=0:numPlots-1
     text(k * xInterval + (xInterval/2) + axisTmp(1) -0.01, axisTmp(3)-1 , H.Label(k+1)) % UNDYNAMISCH !!!!!!!
end;

%XTicks--------------------------------------
xNumTicks = 5;% Fall ohne eventrelated
if ~isempty(strfind(H.type,'eventrelated')) 
    xNumTicks = xNumTicks -1; 
end;

xTicks = [1:((numPlots+1)*xNumTicks)-1]  * xInterval/xNumTicks;  %alle TW durchgehen
xTicks = xTicks + axisTmp(1) ;
set(gca,'XTick',xTicks);

%XTick-Labels

TWlength = 0.150;                       %lnegth of one TWs - HAS TO BE ENTERED MANUALLY
ind = H.T;                               %create a label for each TW
s = size(H.T,2);                         %Number of TWs
ind = ind([1:floor(s/(xNumTicks-1)):s]); %select only as many TWs as xNumTicks
ind = cat(2,ind(2:end) , ind(1));        %place the first element at the end - because labelling starts at the second element
xtickLables = ind - TWlength/2;          %substract TWlength of the center point

xtickLables =  round(xtickLables*10)/10;
set(gca,'XTickLabel',xtickLables);


%YTicks--------------------------------------
ynumTicks = 5;

yTicks = [1:(numPlots+1)*ynumTicks]  * yInterval/ynumTicks; 
yTicks = yTicks +  axisTmp(3); %- yInterval /(2*numYIntervalle)
set(gca,'YTick',yTicks);

%YTick-Labels
%Frequenzzahlen in Array schreiben
ytickLables = ([1:ynumTicks] * ( H.F(end) - H.F(1) ) / ynumTicks ) + H.F(1) ;
for l=1:length(ytickLables)
    ytickLablesStr{l} =  num2str(ytickLables(l)); %Zahlen zu Strings formatieren
end;
ytickLablesStr = repmat(ytickLablesStr,1,numPlots); %den Array so oft wie es Channels gibt aneinanderkopieren
ytickLablesStrLab =  ytickLablesStr;

%an den mittleren TickLabel wird der Name des Chanels vorne ankonkarteniert
for k=0:numPlots-1
    index = k*ynumTicks + floor(ynumTicks/2 +1);
    ytickLablesStrLab{index} = [H.Label{k+1} '    ' ytickLablesStr{k+1}];
end;

set(gca,'YTickLabel',ytickLablesStrLab);

suptitle([ H.method , H.type , ' - Subject ', H.subj , ' - Condition ', H.cond ]);

        end


elseif strcmp(X.datatype,'coupling'),
	clim = [min(0,min(X.data(:))),max(1,max(X.data(:)))];
	ndim = length(size(X.data));
	for k1=1:size(X.data,1)
	for k2=1:size(X.data,2)
		subplot(size(X.data,1),size(X.data,2),(k1-1)*size(X.data,2)+k2);
		if ndim==3,
			area(X.f,squeeze(X.data(k1,k2,:)));
			set(gca,'YLim',clim);
		elseif ndim==4,
			imagesc(X.t,X.f,squeeze(X.data(k1,k2,:,:)));
			set(gca,'CLim',clim);
		end;	
	end;
	end;	

elseif strcmp(X.datatype,'confusion'),
        if nargin>1,
                [kap,sd,H,z,OA,SA,MI]=kappa(X.data);
                fprintf(1,'%s\n',repmat('-',1,8*(size(X.data,1)+1)));
                fprintf(1,'Kappa = %5.3f %c %4.3f(%s)\tOverall Accuracy = %4.1f%%\n',kap,177,sd,repmat('*',sum(-z<norminv([.05,.01,.001]/2)),1),OA*100);
                %disp([X.data,sum(X.data,2);sum(X.data,1),sum(X.data(:))])

                fprintf(1,'%s\n',repmat('-',1,8*(size(X.data,1)+1)));
                for k=1:size(X.data,1),
                        fprintf(1,'%4.1f%%\t',X.data(k,:)./sum(X.data,1)*100);
                        fprintf(1,'| %4.0f\n',sum(X.data(k,:),2));
                end;
                fprintf(1,'%s\n',repmat('-',1,8*(size(X.data,1)+1)));
                fprintf(1,'%4.0f\t',sum(X.data,1));
                fprintf(1,'| %4.0f\n\n',sum(X.data(:)));
        else
                [kap,sd,H,z,OA,SA,MI]=kappa(X.data);
                fprintf(1,'%s\n',repmat('-',1,8*(size(X.data,1)+2)));
                fprintf(1,'Kappa = %5.3f %c %4.3f(%s)\tOverall Accuracy = %4.1f%%\n',kap,177,sd,repmat('*',sum(-z<norminv([.05,.01,.001]/2)),1),OA*100);
                %disp([X.data,sum(X.data,2);sum(X.data,1),sum(X.data(:))])
                fprintf(1,'%s\n',repmat('-',1,8*(size(X.data,1)+2)));

                for k=1:size(X.data,1),
                        fprintf(1,'%6.1f\t',X.data(k,:));
                        fprintf(1,'|%6.1f\t| %4.1f%%\n',sum(X.data(k,:),2),X.data(k,k)/sum(X.data(k,:),2)*100);
                end;
                fprintf(1,'%s\n',repmat('-',1,8*(size(X.data,1)+2)));
                fprintf(1,'%6.0f\t',sum(X.data,1));
                fprintf(1,'|%6.0f\t|\n',sum(X.data(:)));
                fprintf(1,'%s\n',repmat('-',1,8*(size(X.data,1)+1)));
                fprintf(1,'%5.1f%%\t',diag(X.data)'./sum(X.data,1)*100);
                fprintf(1,'|\n\n');

        end;
        if isfield(X,'Label'),
                fprintf(1,'%s\nStage:\t\t',repmat('-',1,8*(size(X.data,1)+2)));
                for k = 1:length(X.Label)
                        fprintf(1,'%-7s\t',X.Label{k});
                end;
        end;
        fprintf(1,'\nSpec.Acc:\t');
        fprintf(1,'%4.1f%%\t',SA*100);
        fprintf(1,'\n%s\n',repmat('-',1,8*(size(X.data,1)+2)));


elseif strcmpi(X.datatype,'pfurtscheller_spectral_difference'),
        nc = ceil(sqrt(X.NS));
        nr = ceil(X.NS/nc);
        nch = size(X.AR,1)/X.NS;
        f = (0:.1:X.SampleRate/2)';
        H = zeros(length(f),X.NC+1);
        for k1=1:nc,
                for k2=1:nr,
                        c = k1+(k2-1)*nc;
                        if nargin>1,
                                H = X.S(:,c+X.NS*(0:X.NC));
                                F = 0:size(X.S,1)-1;
                        else
                                for k3 = 1:X.NC+1;
                                        ix = c + X.NS*(k3-1);
                                        [H(:,k3), F] = freqz(sqrt(X.PE(ix,end)/X.SampleRate),ar2poly(X.AR(ix,:)),f,X.SampleRate);
                                end
                        end;
                        subplot(nc,nr,c);
                        semilogy(F,abs(H),'-');
                        legend({'ref','1','2'});
                        ylabel(sprintf('%s/[%s]^{1/2}',X.PhysDim,X.samplerate_units));
                        v=axis;v(2:4)=[max(F),1e-2,10];axis(v);
                        %hold on;
                        grid on;
                        if isfield(X,'Label');
                                if iscell(X.Label)
                                        title(X.Label{c});
                                else
                                        title(X.Label(c,:));
                                end;
                        else
                                title(['channel # ',int2str(c)]);
                        end;
                end
        end;


elseif strcmpi(X.datatype,'spectrum') || strcmp(X.datatype,'qualitycontrol'),

        if nargin>1,
                Mode=arg2;
        else
                Mode='log';
        end;

       if iscell(X.PhysDim),
              X.PhysDim = strvcat(X.PhysDim);
       end;
       if strcmp(X.datatype,'qualitycontrol'),
                fprintf(1,'\n  [%s]',X.PhysDim(1,:));
                fprintf(1,'\t#%02i',1:size(X.AR,1));
                fprintf(1,'\nMEAN  ');
                fprintf(1,'\t%+7.3f',X.MEAN);
                fprintf(1,'\nRMS');
                fprintf(1,'\t%+7.3f',X.RMS);
                fprintf(1,'\nSTD');
                fprintf(1,'\t%+7.3f',X.STD);
                fprintf(1,'\nQuant');
                fprintf(1,'\t%+7.3f',X.QUANT);
                fprintf(1,'\n  [bit]\nEntropy');
                fprintf(1,'\t%+4.1f',X.ENTROPY);
                fprintf(1,'\n\n');
        end;

        if ~isfield(X,'samplerate_units')
                X.samplerate_units = 'Hz';
        end;
        if ~isfield(X,'PhysDim')
                X.PhysDim = '[1]';
        end;
        if ~isfield(X,'QUANT')
                X.QUANT = 0;
        end;
        if ~isfield(X,'Impedance')
                X.Impedance=5000; %5kOHM
        end;

        if isfield(X,'Label')
                Label = cellstr(X.Label);
        else
                NS = size(X.AR,1);
                Q.Label = [repmat('#',NS,1),int2str([1:NS]')];
        end;

        [n,p] = size(X.AR);
        H=[]; F=[];
        for k=1:size(X.AR,1);
                [h,f] = freqz(sqrt(X.PE(k,size(X.AR,2)+1)/(X.SampleRate*2*pi)),ar2poly(X.AR(k,:)),(0:64*p)/(128*p)*X.SampleRate',X.SampleRate);
                H(:,k)=h(:);F(:,k)=f(:);
        end;
        if strcmp(lower(Mode),'log')
                h=semilogy(F,abs(H),'-',[0,X.SampleRate/2]',[1;1]*mean(X.QUANT)/sqrt(12*X.SampleRate),'k:',[0,X.SampleRate/2]',1e6*[1;1]*sqrt(4*310*138e-25*X.Impedance),'k');
                ylabel(sprintf('%s/[%s]^{1/2}',X.PhysDim(1,:),X.samplerate_units));
                Label = [Label;{'Quantization'};{'Impedance'}];

        elseif strcmp(lower(Mode),'log2')
                semilogy(F,real(H).^2+imag(H).^2,'-',[0,X.SampleRate/2]',[1;1]*X.QUANT.^2/(12*X.SampleRate),'k:');
                ylabel(sprintf('[%s]^2/%s',X.PhysDim(1,:),X.samplerate_units));
                Label = [Label;{'Quantization'}];

        elseif strcmp(lower(Mode),'lin')
                plot(F,abs(H),'-',[0,X.SampleRate/2]',[1;1]*X.QUANT/sqrt(12*X.SampleRate),'k:');
                ylabel(sprintf('%s/[%s]^{1/2}',X.PhysDim(1,:),X.samplerate_units));
                Label = [Label;{'Quantization'}];

        elseif strcmp(lower(Mode),'lin2')
                plot(F,real(H).^2+imag(H).^2,'-',[0,X.SampleRate/2]',[1;1]*X.QUANT.^2/(12*X.SampleRate),'k:');
                ylabel(sprintf('[%s]^2/%s',X.PhysDim(1,:),X.samplerate_units));
                Label = [Label;{'Quantization'}];
        end;
        xlabel(sprintf('f [%s]',X.samplerate_units));
        if isfield(X,'Title'), title(X.Title);
        elseif isfield(X,'FileName'); tmp=X.FileName; tmp(tmp=='_')=' '; title(tmp); end
        H = X;
        legend(Label);

elseif strcmp(X.datatype,'SIESTA_HISTOGRAM')
        if nargin<2,
                chansel=0;
        else
                chansel=arg2;
        end;

        cname=computer;
        if cname(1:2)=='PC',
                PFAD='s:/';
        else
                PFAD='/home/schloegl/';
        end;

        H = load([PFAD,'siesta/t300/',lower(X.filename),'his.mat']);
        R = load([PFAD,'siesta/t300/',lower(X.filename),'res.mat']);

        fn=[PFAD,'siesta/t900/',lower(X.filename),'th.mat'];
        if exist(fn)==2,
                T = load(fn);
        else
                fprintf(2,'Warning: no thresholds available for %s\n',X.filename);
                T = [];
        end;

        H.X = [ones(2^16,1),repmat((-2^15:2^15-1)',1,R.EDF.NS)]*R.EDF.Calib;
        if chansel>0,
                H.H = H.HISTOG(:,chansel);
        else
                H.H = H.HISTOG;
                chansel = 1:R.EDF.NS;
        end;
        H.datatype = 'HISTOGRAM';
        H.N = full(sum(H.H,1));

        if ~isempty(T),
                if any(T.TRESHOLD>=2^15),
                        T.TRESHOLD=T.TRESHOLD-2^15-1;
                        fprintf(2,'Thresholds in %s were by 2^15+1 to high: corrected\n',X.filename);
                end;

                %T.TRESHOLD(:,1)=max(T.TRESHOLD(:,1),R.RES.MU'-7*sqrt(R.RES.SD2'));
                %T.TRESHOLD(:,2)=min(T.TRESHOLD(:,2),R.RES.MU'+7*sqrt(R.RES.SD2'));
        else
                %T.TRESHOLD = ones(R.EDF.NS,1)*[-2^15,2^15-1]; %repmat(nan,R.EDF.NS,2);
                T.TRESHOLD = repmat([2^15-1,-2^15],R.EDF.NS,1)';
                H.Threshold = [ones(2,1),T.TRESHOLD']*R.EDF.Calib(:,chansel);
        end;

        plota(H,'log+');
        suptitle(X.filename);

        
elseif strcmp(X.datatype,'qc:histo')
        if nargin<2,
                yscale = 'log+';
        else
                yscale = arg2;
        end;
        if nargin<3,
                chansel = 0;
        else
                chansel = arg3;
        end;
        if chansel<=0,
                chansel = 1:size(X.HIS.H,2);
        end;

        N = ceil(sqrt(length(chansel)));
        Ny= ceil(length(chansel)/N);
        for k = 1:length(chansel);
        	K = chansel(k); 
                h = X.HIS.H(:,K);
                if isfield(X,'FLAG') && isfield(X.FLAG,'UCAL') && X.FLAG.UCAL, 
       		        t = X.HIS.X(:,min(K,size(X.HIS.X,2)))*X.Calib(K+1,K)+X.Calib(1,K);
        	else
	                t = X.HIS.X(:,min(K,size(X.HIS.X,2)));
                end; 
                if isfield(X,'THRESHOLD'),
	                if isfield(X,'FLAG') && isfield(X.FLAG,'UCAL') && X.FLAG.UCAL, 
	                        MaxMin=X.THRESHOLD(K,[2,1])*X.Calib(K+1,K)+X.Calib(1,K);
	                else        
	                        MaxMin=X.THRESHOLD(K,[2,1]);
			end; 
                elseif isfield(X,'Threshold'),	%%% will become OBSOLETE 
                        MaxMin=X.Threshold(K,:);
                        MaxMin=[max(MaxMin),min(MaxMin)];
                else
                        MaxMin=[max(t) min(t)];
                end;
		h2 = h;
		h2(~xor(t>min(MaxMin),t<max(MaxMin)))=NaN; 

		R.N 	= sumskipnan(h,1);
		R.SUM 	= sumskipnan(h.*t,1);
		R.SSQ 	= sumskipnan(h.*t.*t);
		mu	= R.SUM./R.N;
		R.SSQ0  = R.SSQ-R.SUM.*mu;		% sum square of mean removed
		sd2  	= R.SSQ0./max(R.N-1,0);	     	% variance (unbiased) 
                dT  	= min(abs(diff(t,[],1)));		% QUANT
                xrange  = [min(t(h>0)),max(t(h>0))]; 
                xrange  = xrange + [-1,1]*(diff(xrange)/2+eps);
                
                if strcmp(yscale,'lin '),
                        subplot(Ny,N,k);
                        plot(t,[h],'-');
                elseif strcmp(yscale,'lin+'),
                        subplot(Ny,N,k);
                        plot(t,[h],'-');
                elseif strcmp(yscale,'lin+'),
                        subplot(Ny,N,k);
                        tmp=max(h)/2;
                        tmp=sum(h)/sqrt(2*pi*sd2)*dT/2;
                        %plot(t,[h],'-',t,exp(-(t-mu).^2./sd2/2)./sqrt(2*pi*sd2).*sum(h),'c',mu+sqrt(sd2)*[-5 -3 -1 0 1 3 5],tmp*ones(7,1),'+-',MaxMin,tmp,'rx' );
                        plot(t,[h]+.01,'-',t,exp(-((t-mu).^2)/(sd2*2))/sqrt(2*pi*sd2)*sum(h)*dT,'c',t,h2+.01,'r',mu+sqrt(sd2)*[-5 -3 -1 0 1 3 5],tmp*ones(7,1),'+-',MaxMin,tmp,'rx');
                        v=axis; v=[xrange 1 max(X.HIS.H(:))]; axis(v);
                elseif strcmp(yscale,'stem'),
                        subplot(Ny,N,k);
                        tmp=max(h)/2;
                        tmp=sum(h)/sqrt(2*pi*sd2)*dT/2;
                        stem(t(h>0),h(h>0),'+');
                        hold on
                        plot(t,exp(-((t-mu).^2)/(sd2*2))/sqrt(2*pi*sd2)*sum(h)*dT,'c',t,h2+.01,'r',mu+sqrt(sd2)*[-5 -3 -1 0 1 3 5],tmp*ones(7,1),'+-',MaxMin,tmp,'rx');
                        hold off
                        v=axis; v=[xrange 1 max(X.HIS.H(:))]; axis(v);
                elseif strcmp(yscale,'log ') | strcmp(yscale,'log'),
                        subplot(Ny,N,k);
                        tmp = sqrt(sum(h)/sqrt(2*pi*sd2)*dT);
                        semilogy(t,[h+.01,exp(-((t-mu).^2)/(sd2*2))/sqrt(2*pi*sd2)*sum(h)*dT,h2+.01]);
                elseif strcmp(yscale,'log+'),
                        subplot(Ny,N,k);
                        tmp = sqrt(sum(h(h>0))/sqrt(2*pi*sd2)*dT);
                        semilogy(t,[h+.01,exp(-((t-mu).^2)/(sd2*2))/sqrt(2*pi*sd2)*sum(h(h>0))*dT,h2+.01],'-',mu+sqrt(sd2)*[-5 -3 -1 0 1 3 5]',tmp*ones(7,1),'+-',MaxMin,tmp,'rx');
                        v=axis; v=[xrange 1 max(X.HIS.H(:))]; axis(v);
                elseif strcmp(yscale,'qq'),
                        subplot(Ny,N,k);
                        tmp=.5;sum(h)/2;
                        plot(cumsum(h)/sum(h),normcdf(t,mu,sqrt(sd2)),'xb',[0,1],[0,1],'-c');
                elseif strcmp(yscale,'csum'),
                        subplot(Ny,N,k);
                        tmp=.5;sum(h)/2;
			h2 = cumsum(h); 
			h2((t>min(MaxMin)) & (t<max(MaxMin)))=NaN; 
                        plot(t,cumsum(h)/sum(h),'-',t,normcdf(t,mu,sqrt(sd2)),'c',t,h2,'r',mu+sqrt(sd2)*[-5 -3 -1 0 1 3 5]',tmp*ones(7,1),'+-',MaxMin,tmp,'rx');
                        v=axis; v(1:2)=[MaxMin(2)+0.1*diff(MaxMin) MaxMin(1)-0.1*diff(MaxMin)]; axis(v);
                elseif strcmp(yscale,'CDF'),
                        subplot(Ny,N,k);
                        tmp = sum(h)/2;
			cdf = cumsum(h)/sum(h); 
			h2  = cdf; 
			h2((t>min(MaxMin)) & (t<max(MaxMin)))=NaN; 
			plot(t,cdf,'-b',t,h2,'-r',mu+sqrt(sd2)*[-5 -3 -1 0 1 3 5]',ones(1,7)/2,'-+g',MaxMin,.5,'rx',[min(t),max(t)],[0,1],'bx'); 
			grid on; 
                        %plot([X.X(1,:)-eps;X.X],[zeros(1,size(X.H,2));cumsum(X.H,1)]./X.N(ones(size(X.X,1)+1,1),:),'-');
                        t = [t(1)-eps;t];
                        v = axis; v(3:4) = [0,1]; axis(v);
                elseif strcmp(yscale,'stacked'),
                        bar(t,h,'stacked');
                end;
                if isfield(X,'Label');
	                title(X.Label{K});
                end; 
                if isfield(X,'PhysDim');
                	xlabel(['[',X.PhysDim{K},']']);
                end; 	
        end;
        if iscell(X.PhysDim),
%               X.PhysDim = strvcat(X.PhysDim);
        end;
       	X.RES = hist2res(X); %this uses the scaling 

                fprintf(1,'\nLabel:');
                fprintf(1,'\t%7s',X.Label{:});
                fprintf(1,'\nUnits:');
                fprintf(1,'\t[%5s]',X.PhysDim{:});
                fprintf(1,'\nMEAN:');
                fprintf(1,'\t%+7.3f',X.RES.MEAN);
                fprintf(1,'\nRMS:');
                fprintf(1,'\t%+7.3f',X.RES.RMS);
                fprintf(1,'\nSTD:');
                fprintf(1,'\t%+7.3f',X.RES.STD);
                fprintf(1,'\nQuant:');
                fprintf(1,'\t%+7.3f',X.RES.QUANT);
                fprintf(1,'\nEntropy [bit]:\n');
                fprintf(1,'\t%+4.1f',X.RES.ENTROPY);
                fprintf(1,'\nRatio of lost samples [%%]:\n');
                fprintf(1,'\t%6.3f',X.RES.ratio_lost*100);
                fprintf(1,'\n\n');

        if ~isfield(X,'AR')
        	return; 
       	end; 
        fprintf(1,'\t#%02i',1:size(X.AR,1));
        [n,p] = size(X.AR);
        H=[]; F=[];
        for k=1:size(X.AR,1);
                [h,f] = freqz(sqrt(X.PE(k,size(X.AR,2)+1)/(X.SampleRate*2*pi)),ar2poly(X.AR(k,:)),(0:64*p)/(128*p)*X.SampleRate',X.SampleRate);
                H(:,k)=h(:);F(:,k)=f(:);
        end;
        Mode ='log';
        if ~isfield(X,'Impedance')
        	X.Impedance = 5000; 
        	warning('Impedance not available assume 5kOhm');
        end; 
        if ~isfield(X,'samplerate_units')
		X.samplerate_units = 'Hz';
	end;
	if ~isfield(X,'Label')
		Label = [repmat('ch ',n,1);int2str([1:n]')];
	else
		Label = X.Label;	
	end;
	
	figure(2);
        if strcmp(lower(Mode),'log')
                h=semilogy(F,abs(H),'-',[0,X.SampleRate/2]',[1;1]*mean(X.RES.QUANT)/sqrt(12*X.SampleRate),'k:',[0,X.SampleRate/2]',1e6*[1;1]*sqrt(4*310*138e-25*X.Impedance),'k');
                ylabel(sprintf('%s/[%s]^{1/2}',X.PhysDim(1,:),X.samplerate_units));
                Label = [Label;{'Quantization'};{'Impedance'}];

        elseif strcmp(lower(Mode),'log2')
                semilogy(F,real(H).^2+imag(H).^2,'-',[0,X.SampleRate/2]',[1;1]*X.RES.QUANT.^2/(12*X.SampleRate),'k:');
                ylabel(sprintf('[%s]^2/%s',X.PhysDim(1,:),X.samplerate_units));
                Label = [Label;{'Quantization'}];

        elseif strcmp(lower(Mode),'lin')
                plot(F,abs(H),'-',[0,X.SampleRate/2]',[1;1]*X.RES.QUANT/sqrt(12*X.SampleRate),'k:');
                ylabel(sprintf('%s/[%s]^{1/2}',X.PhysDim(1,:),X.samplerate_units));
                Label = [Label;{'Quantization'}];

        elseif strcmp(lower(Mode),'lin2')
                plot(F,real(H).^2+imag(H).^2,'-',[0,X.SampleRate/2]',[1;1]*X.RES.QUANT.^2/(12*X.SampleRate),'k:');
                ylabel(sprintf('[%s]^2/%s',X.PhysDim(1,:),X.samplerate_units));
                Label = [Label;{'Quantization'}];
        end;
        xlabel(sprintf('f [%s]',X.samplerate_units));
        if isfield(X,'Title'), title(X.Title);
        elseif isfield(X,'FileName'); tmp=X.FileName; tmp(tmp=='_')=' '; title(tmp); end
        H = X;
        legend(Label);
        
        
elseif strcmp(X.datatype,'HISTOGRAM') || strcmp(X.datatype,'qc:histo')
        
        if nargin<3,
                chansel=0;
        else
                chansel=arg3;
        end;
        if nargin<2
                yscale='lin ';
        else
                yscale=arg2;
        end;

        if ~isfield(X,'N');
                X.N = sumskipnan(full(X.H),1);
        end;

        if chansel<=0,
                chansel = 1:size(X.H,2);
        end;

        N=ceil(sqrt(size(X.H,2)));
        for K = chansel;
                %min(K,size(X.X,2))
                t = X.X(:,min(K,size(X.X,2)));
                %HISTO=hist2pdf(HISTO);
                h = X.H(:,K);

                mu = (t(h>0)'*h(h>0))/X.N(K);%sumskipnan(repmat(t,size(h)./size(t)).*h,1)./sumskipnan(h,1);
                x  = t-mu; %(repmat(t,size(h)./size(t))-repmat(mu,size(h)./size(mu)));
                sd2= sumskipnan((x(h>0).^2).*h(h>0),1)./X.N(K);

                [tmp,tmp2]=find(h>0);

                if isfield(X,'Threshold'),
                        MaxMin=X.Threshold(:,K)';
                        MaxMin=[max(MaxMin),min(MaxMin)];
                else
                        MaxMin=t([max(tmp) min(tmp)]);
                end;

                if strcmp(yscale,'lin '),
                        subplot(ceil(size(X.H,2)/N),N,K);
                        plot(t,[h],'-');
                elseif strcmp(yscale,'lin+'),
                        subplot(ceil(size(X.H,2)/N),N,K);
                        tmp = diff(t);
                        dQ  = 1;min(tmp(tmp>0));
                        tmp=max(h)/2;
                        tmp=sum(h)/sqrt(2*pi*sd2)*dQ/2;
                        %plot(t,[h],'-',t,exp(-(t-mu).^2./sd2/2)./sqrt(2*pi*sd2).*sum(h),'c',mu+sqrt(sd2)*[-5 -3 -1 0 1 3 5],tmp*ones(7,1),'+-',MaxMin,tmp,'rx' );
                        plot(t,[h]+.01,'-',t,exp(-((t-mu).^2)/(sd2*2))/sqrt(2*pi*sd2)*sum(h)*dQ,'c',mu+sqrt(sd2)*[-5 -3 -1 0 1 3 5],tmp*ones(7,1),'+-',MaxMin,tmp,'rx');
                        v=axis; v=[MaxMin(2) MaxMin(1) 1 max(h)]; axis(v);
                elseif strcmp(yscale,'log ') | strcmp(yscale,'log'),
                        subplot(ceil(size(X.H,2)/N),N,K);
                        tmp = diff(t);
                        dQ  = min(tmp(tmp>0));
                        tmp = sqrt(sum(h)/sqrt(2*pi*sd2)*dQ);
                        %semilogy(t,[h],'-')
                        semilogy(t,[h+.01,exp(-((t-mu).^2)/(sd2*2))/sqrt(2*pi*sd2)*sum(h)*dQ]);
                elseif strcmp(yscale,'log+'),
                        subplot(ceil(size(X.H,2)/N),N,K);
                        tmp = diff(t);
                        if any(tmp>0)
	                        dQ = min(tmp(tmp>0));
	                else 
	                	dQ = 1;
	                end;	        
                        tmp = sqrt(sum(h(h>0))/sqrt(2*pi*sd2)*dQ);
                        %semilogy(t,[h]+.01,'-',t,exp(-(t*ones(size(mu))-ones(size(t))*mu).^2./(ones(size(t))*sd2)/2)./(ones(size(t))*(sqrt(2*pi*sd2)./sum(h))),'c',mu+sqrt(sd2)*[-5 -3 -1 0 1 3 5],tmp*ones(7,1),'+-',MaxMin,tmp,'rx');
                        %semilogy(t,[h]+.01,'-',t,exp(-(t(:,ones(size(mu)))-mu(ones(size(t)),:)).^2./sd2(ones(size(t)),:)/2)./sqrt(2*pi*sd2(ones(size(t)),:)).*(ones(size(t))*sum(h)),'c',mu+sqrt(sd2)*[-5 -3 -1 0 1 3 5],tmp*ones(7,1),'+-',MaxMin,tmp,'rx');
                        %semilogy(t,[h]+.01,'-',t,exp(-((t-mu).^2)/(sd2*2))/sqrt(2*pi*sd2)*sum(h)*dQ,'c',mu+sqrt(sd2)*[-5 -3 -1 0 1 3 5]',tmp*ones(7,1),'+-',MaxMin,tmp,'rx');
                        semilogy(t,[h+.01,exp(-((t-mu).^2)/(sd2*2))/sqrt(2*pi*sd2)*sum(h(h>0))*dQ],'-',mu+sqrt(sd2)*[-5 -3 -1 0 1 3 5]',tmp*ones(7,1),'+-',MaxMin,tmp,'rx');
                        v=axis; v=[MaxMin(2)+0.1*diff(MaxMin)-50*eps MaxMin(1)-0.1*diff(MaxMin)+50*eps 1 max(h)];
                        axis(v);
                        v=axis; v=[v(1:2) 1 max(h)]; axis(v);
                elseif strcmp(yscale,'qq'),
                        subplot(ceil(size(X.H,2)/N),N,K);
                        tmp=.5;sum(h)/2;
                        %plot(t,cumsum(h)/sum(h),'-',t,cumsum(exp(-(t-mu).^2/sd2/2)/sqrt(2*pi*sd2)/X.N(K)),'c',mu+sqrt(sd2)*[-5 -3 -1 0 1 3 5]',tmp*ones(7,1),'+-',MaxMin,tmp,'rx');
                        plot(cumsum(h)/sum(h),normcdf(t,mu,sqrt(sd2)),'xb',[0,1],[0,1],'-c');
                elseif strcmp(yscale,'csum'),
                        subplot(ceil(size(X.H,2)/N),N,K);
                        tmp=.5;sum(h)/2;
                        %plot(t,cumsum(h)/sum(h),'-',t,cumsum(exp(-(t-mu).^2/sd2/2)/sqrt(2*pi*sd2)/X.N(K)),'c',mu+sqrt(sd2)*[-5 -3 -1 0 1 3 5]',tmp*ones(7,1),'+-',MaxMin,tmp,'rx');
                        plot(t,cumsum(h)/sum(h),'-',t,normcdf(t,mu,sqrt(sd2)),'c',mu+sqrt(sd2)*[-5 -3 -1 0 1 3 5]',tmp*ones(7,1),'+-',MaxMin,tmp,'rx');
                        v=axis; v(1:2)=[MaxMin(2)+0.1*diff(MaxMin) MaxMin(1)-0.1*diff(MaxMin)]; axis(v);
                elseif strcmp(yscale,'CDF'),
                        %subplot(ceil(size(X.H,2)/N),N,K);
                        tmp=sum(h)/2;
                        %semilogx(X.X,cumsum(X.H,1)./X.N(ones(size(X.X,1),1),:),'-');
                        plot([X.X(1,:)-eps;X.X],[zeros(1,size(X.H,2));cumsum(X.H,1)]./X.N(ones(size(X.X,1)+1,1),:),'-');
                        t = [t(1)-eps;t];
                        %plot(t,[cumsum([0;h])/X.N(K),normcdf(t,mu,sqrt(sd2))])
                        %plot(t,cumsum([0;h])/X.N(K))
                        %v=axis; v(1:2)=[MaxMin(2)+0.1*diff(MaxMin) MaxMin(1)-0.1*diff(MaxMin)]; axis(v);
                        v=axis; v(3:4)=[0,1]; axis(v);
                elseif strcmp(yscale,'stacked'),
                        bar(t,h,'stacked');
                end;
        end;
elseif strcmp(X.datatype,'DBI-EPs'),
        if nargin<2,
                arg2='b';
        end;
        if nargin<3,
                arg3=100;
        end;

        for k=1:length(X.eegchan),
                subplot(13,10,k);

                mu=X.RES(k).SUM./X.RES(k).N;
                sd=sqrt((X.RES(k).SSQ-X.RES(k).SUM.*mu)./max(X.RES(k).N,0));
                se=sd./sqrt(X.RES(k).N);

                h=plot(X.t,[mu(:),sd(:),se(:)]*[1,1,1,1,1;-1,1,0,0,0;0,0,-1,1,0],arg2);
                set(h(5),'linewidQh',2);
                set(h(1),'color',[0,.75,.75]);
                set(h(2),'color',[0,.75,.75]);
                set(h(3),'color',[.5,.5,1]);
                set(h(4),'color',[.5,.5,1]);
                axis([-3,3,-arg3,arg3]);
                title(sprintf('#%i',X.eegchan(k)));
        end;


elseif strcmp(X.datatype,'SCATTER'),
        s='x';

        if nargin<2,

        elseif nargin==2,
                if length(arg2)<3,
                        s = arg2;
                else
                        Labels=arg2;
                end
        elseif nargin>2,
                s = arg2;
                Labels = arg3;
        end;
	if ~exist('Labels','var') 
        	Labels = cellstr(int2str([1:size(X.data,2)]'));
        end; 	
        if isfield(X,'Classlabel')
        	CL = unique(X.Classlabel);
        else 
        	CL = 1; 	
        end;	

        if length(X)==1,
                if ~isfield(X,'R');
                        [X.R,X.p,X.CIL,X.CIU] = corrcoef(X.data,'Rank');
                end;
                [nr,nc] = size(X.data);
                nc=nc-1;
                for k   = 1:nc,
                        for l = k+1:nc+1,%[1:k-1,k+1:nc],
                                h=subplot(nc,nc,(k-1)*nc+l-1);
                                if length(CL)==1,
	                                plot(X.data(:,l),X.data(:,k),s);
                                elseif length(CL)==2,
	                                plot(X.data(X.Classlabel==CL(1),l),X.data(X.Classlabel==CL(1),k),'bx',X.data(X.Classlabel==CL(2),l),X.data(X.Classlabel==CL(2),k),'ro');
	                        end;        
                                ht=title(sprintf('r=%.3f %s [%.3f,%.3f]',X.R(k,l),char('*'*(X.p(k,l)<[.05,.01,.001])),X.CIL(k,l),X.CIU(k,l)));
                                %ht=title(sprintf('r = %.4f %s',X.R(k,l),char('*'*(X.p(k,l)<[.05,.01,.001]))));
                                pos=get(ht,'position');
                                %pos(2)=max(X.data(k));
                                %set(ht,'position',pos);
                                % title gives the r-value, its confidence interval (see CORRCOFF for which level),
                                %      and indicates the significance for alpha=0.05 (*), 0.01(**) and 0.001(***)
                                if l == (k+1),
                                        xlabel(Labels{l});
                                        ylabel(Labels{k});
                                else
                                        set(h,'xtick',[]);
                                        set(h,'ytick',[]);
                                end;
                        end;
                end;
        else
                if ~isfield(X,'R');
                        [X.R,X.p,X.CIL,X.CIU] = corrcoef(X(1).data,X(2).data,'Rank');
                end;
                [nr,nc2] = size(X(2).data);
                for k   = 1:nc,
                        for l = 1:nc2,%[1:k-1,k+1:nc],
                                h=subplot(nc,nc2,(k-1)*nc2+l-1);
                                plot(X(2).data(:,l),X(1).data(:,k),s);
                                ht=title(sprintf('r=%.3f %s [%.3f,%.3f]',X.R(k,l),char('*'*(X.p(k,l)<[.05,.01,.001])),X.CIL(k,l),X.CIU(k,l)));
                                %ht=title(sprintf('r = %.4f %s',X.R(k,l),char('*'*(X.p(k,l)<[.05,.01,.001]))));
                                pos=get(ht,'position');
                                %pos(2)=max(X.data(k));
                                %set(ht,'position',pos);
                                % title gives the r-value, its confidence interval (see CORRCOFF for which level),
                                %      and indicates the significance for alpha=0.05 (*), 0.01(**) and 0.001(***)
                                if l == 1,
                                        %xlabel(arg3{l});
                                        ylabel(arg3{k});
                                else
                                        %set(h,'xtick',[]);
                                        set(h,'ytick',[]);
                                end;
                        end;
                        if k == nc,
                                xlabel(arg3{l});
                                %ylabel(arg3{k});
                        else
                                set(h,'xtick',[]);
                                %set(h,'ytick',[]);
                        end;
                end;
        end;

elseif strcmp(X.datatype,'STAT2'),
        if nargin<2,
                arg2='b';
        end;
        if isfield(X,'t')
                t=X.t;
        else
                t=1:length(X.SUM);
        end;
        if nargin>2
                t=arg3;
        end;

        mu=X.SUM./X.N;
        sd=sqrt((X.SSQ-X.SUM.*mu)./max(X.N-1,0));
        se=sd./sqrt(X.N);

        h=plot(t,[mu(:),sd(:),se(:)]*[1,1,1,1,1;0,-1,1,0,0;0,0,0,-1,1],arg2);
        set(h(1),'linewidQh',2);
        tmp=get(h(1),'color');
        set(h(2),'color',1-(1-tmp)/2);
        set(h(3),'color',1-(1-tmp)/2);
        %set(h(4),'color',[0,0,1]);
        %set(h(5),'color',[0,0,1]);


elseif strcmp(X.datatype,'TSD1'),
        if nargin<2,
                arg2='b';
        end;
        h=plot(X.t,[X.mu,X.sd]*[1,1,1;0,-1,1],arg2);
        set(h(1),'linewidQh',2);
        hold on
        h=plot(X.TI(:),1,'.k');


elseif strcmp(X.datatype,'MEAN+STD')

        if isfield(X,'MEAN');
                sz = [size(X.MEAN),1];
        elseif isfield(X,'SUM');
                sz = [size(X.SUM),1];
        end;
        nchns = sz(2);  % Number of channels

        if nargin < 2
                clf;
                nf = [];
        else
                nf = arg2;  % Handles to subplots
        end;
        sz,
        if isempty(nf)
                for k0 = 1:sz(3),
                        if sz(3)>1, figure(k0); end;
                        for k = 1:nchns
                                nf(k0,k) = subplot(ceil(nchns/ceil(sqrt(nchns))),ceil(sqrt(nchns)),k);  % Handles to subplots
                        end;
                end;
        end;

        if isfield(X,'Label')
                if ischar(X.Label)
                        X.Label = cellstr(X.Label);
                end;
        end;

        if (~isfield(X,'MEAN') || ~isfield(X,'STD')) %& isfield(X,'SUM') & isfield(X,'N') & isfield(X,'SSQ')
                X.MEAN 	= X.SUM./X.N;			% mean
                X.MSQ  	= X.SSQ./X.N;;			% mean square
                X.RMS  	= sqrt(X.MSQ);			% root mean square
                %X.SSQ0	= X.SSQ-X.SUM.*X.MEAN;		% sum square of mean removed
                X.SSQ0	= X.SSQ - real(X.SUM).*real(X.MEAN) - imag(X.SUM).*imag(X.MEAN);	% sum square of mean removed
                n1 	= max(X.N-1,0);			% in case of n=0 and n=1, the (biased) variance, STD and SEM are INF
                X.VAR  	= X.SSQ0./n1;	     		% variance (unbiased)
                X.STD  	= sqrt(X.VAR);		     	% standard deviation
                X.SEM  	= sqrt(X.SSQ0./(X.N.*n1)); 	% standard error of the mean
        end;

        if nargin>2
                minmean = arg3;
                maxmean = arg4;
                maxstd  = arg5;
        else
                minmean = floor(min(X.MEAN(:)));
                maxmean = ceil(max(X.MEAN(:)));
                maxstd  = ceil(max(X.STD(:)));
        end;

        for k0 = 1:sz(3),
                if sz(3)>1, figure(k0); end;
                for k = 1:sz(2)  % For each channel

                        subplot(nf(k0,k));
                        [ax,h1,h2] = plotyy(X.T,X.MEAN(:,k,k0),X.T,X.STD(:,k,k0));
                        drawnow;
                        set(ax(1),'FontSize',8);
                        set(ax(2),'FontSize',8);

                        % Sets the axes limits to avoid overlapping of the two functions
                        set(ax(1),'YLim',[minmean-maxstd maxmean]);
                        set(ax(2),'YLim',[0 -minmean+maxstd+maxmean]);
                        set(ax,'XLim',[min(X.T) max(X.T)]);

                        set(ax,'YTickLabel',[]);
                        set(ax,'YTick',[]);

                        % Label y1-axis (mean)
                        temp = [floor(minmean/10) * 10 : 10 : ceil(maxmean/10) * 10];  % Label only ..., -20, -10, 0, 10, 20, 30, ...
                        temp = [0:5]/5*(maxmean-minmean)+minmean; % [floor(minmean/10) * 10 : 10 : ceil(maxmean/10) * 10];  % Label only ..., -20, -10, 0, 10, 20, 30, ...
                        set(ax(1),'YTick',temp);

                        set(ax(1),'YTickLabel',temp);

                        % Label y2-axis (standard deviation)
                        temp = [0 : 5] * ceil(maxstd/5);  % Label only 0, 10, 20, 30, ...
                        set(ax(2),'YTick',temp);

                        set(ax(2),'YTickLabel',temp);

                        % Label x-axis

                        xlabel('Time (s)');
                        grid on;

                        if isfield(X,'Label')  % Print label of each channel (if such a label exists)
                                if k <= length(X.Label)
                                        title(X.Label{k},'FontSize',8,'Interpreter','none');
                                end;
                        end;

                        if isfield(X,'trigger')  % Mark trigger
                                line([X.T(X.trigger),X.T(X.trigger)],[minmean-maxstd,maxmean],'Color',[1 0 0]);
                        end;
                end;
        end;
        drawnow;
        h = ax;
        %set(0,'DefaultTextInterpreter','none');  % Avoid having TeX interpretation in title string
        %suptitle(X.Title);


elseif strcmp(X.datatype,'CORRELATION_WITH_REFERENCE')
        if isfield(X,'ELEC') | isfield(X,'ELPOS')
                fprintf(2,'PLOTA: X.ELPOS, X.ELEC not supported yet.\n');
        end;
        FLAG.TOPOMAP = 1;
        if nargin>1,
                if ~exist(arg2,'file');
                        fprintf(2,'Warning PLOTA: electrode position file not found.\n');
                        FLAG.TOPOMAP = 0;
                end;
        else
                %	fprintf(2,'Warning PLOTA: electrode position file not specified.\n');
                FLAG.TOPOMAP = 0;
        end;
        if FLAG.TOPOMAP,
                topoplot(X.corr,arg2,'maplimits',[-1,1]);
                colorbar;
        else
                plot(X.corr);
                axis([1,length(X.corr),-1,1]);
        end;


elseif strcmp(X.datatype,'Classifier')
        if ~isfield(X,'tsc'),
                X.tsc=X.TI*16+[-15,1];
        end;
        if (nargin==1);
                arg2 = 'all';
                for k=1:4,
                        hf(k) = subplot(1,4,k);
                end
        elseif (nargin==2);
                if isnumeric(arg2)
                        hf = arg2;
                        if length(hf)==3,
                                arg2='all';
                        elseif length(hf)==1,
                                arg2='acc';
                                subplot(hf);
                        end;
                else   % arg2=arg2;
                end
        elseif (nargin==3);
                if isnumeric(arg2)
                        hf = arg2;
                        arg2 = arg3;
                else
                        hf = arg3;
                end;
        else

        end;

        if ~isfield(X,'T');
                %X.T = (1:size(X.acc,1))';
                X.T = (1:size(X.MDA.ACC00,1))';
                if isfield(X,'Fs'),
                        X.T = X.T/X.Fs;
                end;
        else;

        end;
        LEG = [];
        if isfield(X,'Classes'),
                if ischar(X.Classes)
                        LEG = X.Classes;
                elseif isnumeric(X.Classes)
                        LEG = num2str(X.Classes(:));
                end;
        end;
        if strncmpi(arg2,'acc',3)
                if isfield(X,'tsc'),
                        patch(X.T(X.tsc([1,1,2,2,1])),[0,1,1,0,0]*100,[1,1,1]*.8);
                end;
                hold on;
                %plot(X.T,X.acc*100,'-',X.T([1,end]),[100,100]./size(X.acc,2),'k:');
                plot(X.T,X.MDA.ACC00*100,'-',X.T([1,end]),[100,100]./size(X.MDA.ACC00,2),'k:');
                hold off;
                v=axis;v(3:4)=[0,100];axis(v);

                ylabel('Accuracy [%]');
                grid on;
                if ~isempty(LEG)
                        legend(LEG);
                end

        elseif strcmpi(arg2,'fixed') |  strcmpi(arg2,'fixed-MDA'),
                if 0,isfield(X,'tsc'),
                        patch(X.T(X.tsc([1,1,2,2,1])),[0,1,1,0,0]*100,[1,1,1]*.8);
                end;
                hold on;
                plot(X.T,X.MDA.acc*100,'-',X.T([1,end]),[100,100]./size(X.MDA.acc,2),'k:');
                hold off;
                v=axis;v(3:4)=[0,100];axis(v);

                ylabel('mean recognistion rate [%]');
                grid on;
                if ~isempty(LEG)
                        legend(LEG);
                end

        elseif strcmpi(arg2,'fixed-LLH')
                if 0,isfield(X,'tsc'),
                        patch(X.T(X.tsc([1,1,2,2,1])),[0,1,1,0,0]*100,[1,1,1]*.8);
                end;
                hold on;
                plot(X.T,X.LLH.acc*100,'-',X.T([1,end]),[100,100]./size(X.LLH.acc,2),'k:');
                hold off;
                v=axis;v(3:4)=[0,100];axis(v);

                ylabel('Accuracy [%]');
                grid on;
                if ~isempty(LEG)
                        legend(LEG);
                end

        elseif strcmpi(arg2,'KAPPA') | strcmpi(arg2,'KAPPA-MDA')
                if isfield(X,'tsc'),
                        patch(X.T(X.tsc([1,1,2,2,1])),[0,1,1,0,0]*100,[1,1,1]*.8);
                end;
                hold on;
                plot(X.T,[X.MDA.KAP00,X.MDA.ACC00]*100);
                hold off;
                grid on;
                %v=axis; v(3:4)=[-10,100]; axis(v);
                v=axis; v(3:4)=[0,100]; axis(v);
                ylabel('Kappa [%], Accuracy [%]')
                xlabel('time t [s]');
                legend('Kappa', 'Accuracy');

        elseif strcmpi(arg2,'KAPPA-LLH')
                if isfield(X,'tsc'),
                        patch(X.T(X.tsc([1,1,2,2,1])),[0,1,1,0,0]*100,[1,1,1]*.8);
                end;
                hold on;
                plot(X.T,[X.LLH.KAP00,X.LLH.ACC00]*100);
                hold off;
                grid on;
                %v=axis; v(3:4)=[-10,100]; axis(v);
                v=axis; v(3:4)=[0,100]; axis(v);
                ylabel('Kappa [%], Accuracy [%]')
                xlabel('time t [s]');
                legend('Kappa', 'Accuracy');

        elseif strncmpi(arg2,'MI',2)
                if isfield(X,'tsc'),
                        patch(X.T(X.tsc([1,1,2,2,1])),[0,1,1,0,0],[1,1,1]*.8);
                end;
                hold on;
                if isfield(X,'I0') %& any(X.I0(:)~=0);
                        h=plot(X.T,[X.I0,X.I])
                else
                        h=plot(X.T,[X.MD2.I0,sum(X.MD2.I0,2)]);
                        %h=plot(X.T,[X.GRB.I0,X.GRB.I2]);
                end;
                hold off
                grid on;
                v=axis; v(3:4)=[0,1]; axis(v);
                set(h(end),'linewidQh',2)
                ylabel('MI [bit]');
                xlabel('time t [s]');
                if ~isempty(LEG)
                        legend(LEG);
                end

        elseif strncmpi(arg2,'all',3)
                plota(X,'acc',hf(1));
                plota(X,'KAPPA',hf(2));
                plota(X,'fixed',hf(3));
                plota(X,'MI',hf(4));
        end;

elseif strncmp(X.datatype,'TSD_BCI',7) && (nargin>1) && strcmpi(arg2,'TSD');
        N = length(X.CL);
        if N>2,
                for k=1:N,
                        nf(k)=subplot(ceil(sqrt(N)),N/ceil(sqrt(N)),k);
                end;
        end;

        if N==2, N=1; end;
        for k=1:N,
                if N>1,
                        subplot(nf(k));
                end;
                h=plot(X.T,[X.MEAN1(:,k),X.MEAN2(:,k),X.SD1(:,k),X.SD2(:,k)]*[1,0,1,0,1,0; 0,1,0,1,0,1; 0,0,1,0,-1,0; 0,0,0,1,0,-1]);
                set(h(1),'linewidQh',2);
                set(h(2),'linewidQh',2);
                mset(h(3:6),'linewidQh',1);
                mset(h([1,3,5]),'color','b');
                mset(h([2,4,6]),'color','g');
                ylabel('Average TSD');
                v=axis;v(1:2)=[min(X.T),max(X.T)];axis(v);
        end;


elseif isfield(X,'TSD') && isfield(X.TSD,'datatype') && strcmp(X.TSD.datatype,'TSD_BCI9')
        if nargin<2,
                clf;
                for k=1:6,
                        nf(k)=subplot(3,2,k);
                end;
        else
                nf=arg2;
        end;

        fid = 1;
       	tix = X.TSD.tix(1);
        N = length(X.TSD.CL);
       	c = (N-1)/N;

	fprintf(fid,'Classifier:	%s\n',X.datatype);
	if isfield(X,'hyperparameters')
		f = fieldnames(X.hyperparameters); 
		for k=1:length(f),
			fprintf(fid,'	%s = %f\n',f{k},getfield(X.hyperparameters,f{k}));
		end; 	
	end; 	
        fprintf(fid,'Error:   		%4.1f %% \n',100-X.TSD.ACC00(tix)*100);
       	fprintf(fid,'Accuracy:		%4.1f %% \nspecific Accuracy:	',X.TSD.ACC00(tix)*100);
       	[kap,sd,H,z,OA,SA,MI] = kappa(X.TSD.optCMX);
       	fprintf(fid,'%4.1f   ',SA*100);
       	fprintf(fid,'\nKappa:   		%4.2f %c %4.2f\n',X.TSD.KAP00(tix),177,X.TSD.Ksd00(tix));
       	fprintf(fid,'I(Wolpaw):		%4.2f bit\n',wolpaw_entropy(X.TSD.ACC00(tix),N));
       	fprintf(fid,'I(Nykopp):		%4.2f bit\n',X.TSD.I_Nykopp(tix));

        if isfield(X.TSD,'I'),
	        fprintf(fid,'I(Continous):     SUM = %4.2f  [ ',sumskipnan(X.TSD.I(tix,:))*c);
        	fprintf(fid,'%4.2f   ',X.TSD.I(tix,:));
        	t = X.T.t-X.T.t0; t(t<.5)=NaN;
        	fprintf(fid,' ] \nSTMI:   		%4.2f  [ ',max([sum(X.TSD.I,2)]./t)*c);
        	fprintf(fid,'   %4.2f',max(X.TSD.I./[t(:,ones(1,size(X.TSD.I,2)))]));
        	fprintf(fid,' ] \nSNR:			')
        	fprintf(fid,'%4.2f   ',X.TSD.SNR(tix,:));
        	fprintf(fid,'\ncorrelation (parametric):	')
        	fprintf(fid,'%4.2f   ',X.TSD.r(tix,:));
        	fprintf(fid,'\nrank correlation:	  ');
	        fprintf(fid,'\nAUC:			');
	        fprintf(fid,'%4.2f   ',X.TSD.AUC(tix,:));
	        fprintf(fid,'\n');
        end; 
        if isfield(X,'rankcorrelation');
                fprintf(fid,'%4.2f   ',X.TSD.rankcorrelation(tix,:));
        end;

        xlim = [min(0,X.T.t(1)),max(X.T.t(end))];
        
        subplot(nf(1));
        if ~isfield(X.TSD,'Labels')
                for k = 1:length(X.TSD.CL), Labels{k}=int2str(X.TSD.CL(k)); end;
        else
                Labels = X.TSD.Labels;
        end;
        %plota(X);

        plot(X.T.t, 100*[X.TSD.ACC00, X.TSD.KAP00])
        axis([xlim,-20,100])
        xlabel('time [s]')
        title('Accuracy and Kappa')
        legend({'Accuracy [%]','Kappa [%]'})

        subplot(nf(3));
        plot(X.T.t,[sum(X.TSD.I,2)*c,X.TSD.I_Nykopp(:),wolpaw_entropy(X.TSD.ACC00,N)])
        legend({'I_{continous}','I_{Nykopp}','I_{Wolpaw}'})
        title('Mutual information')
        ylabel('I [bit]');
        xlabel('time [s]')
        v=axis; axis([xlim,0,1.5]);

        if 0,
                subplot(nf(5));
                plot(X.T.t,X.MEAN2)
                hold on
                plot(X.T.t,X.MEAN1)
                hold off
                Title('average output of one-vs-rest classifiers')
                ylabel('TSD')
                xlabel('time [s]')
                tmp = strvcat(Labels);
                legend(cellstr([tmp,repmat('(+)',size(tmp,1),1);tmp,repmat('(-)',size(tmp,1),1)]))
                v=axis; axis([xlim,v(3:4)]);
        else
                subplot(nf(5));
                t = X.T.t-X.T.t0;
                %t(t < 3.5)=NaN;
                t(t < 0.5)=NaN;
                plot(X.T.t,[sum(X.TSD.I,2)*c,X.TSD.I_Nykopp(:),wolpaw_entropy(X.TSD.ACC00,N)]./repmat(t,1,3))
                legend({'STMI_{C}','STMI_{N}','STMI_{W}'})
                title('Steepness of Mutual information')
                ylabel('STMI [bit/s]');
                xlabel('time [s]')
                v=axis; axis([xlim,0,1]);
        end;

        subplot(nf(2));
        plot(X.T.t,X.TSD.AUC)
        xlabel('time [s]')
        ylabel('AUC [1]')
        title('area-under-the-(ROC)-curve');
        legend(Labels)
        v=axis; axis([xlim,.4,1]);

        subplot(nf(4));
        plot(X.T.t,[X.TSD.I,sum(X.TSD.I,2)*c,])
        xlabel('time [s]')
        ylabel('I [bit]');
        title('Mutual Information of continous output');
        legend([Labels,{'SUM'}])
        v=axis; axis([xlim,0,1.5]);

        subplot(nf(6));
        if isfield(X.TSD,'N')
        	% show significance interval 
	        alpha = .05;
 	       	ci = tanh(sqrt(2)*erfinv(1-2*alpha)./sqrt(X.TSD.N-3));		% confidence interval for alpha of z
	        plot(X.T.t,X.TSD.r,'-',X.T.t,ci*[-1,1],'k:');
	        LEG = [Labels,{sprintf('alpha=%f',alpha)}];
	else
	        plot(X.T.t,X.TSD.r,'-');
	        LEG = Labels;
        end; 
        xlabel('time [s]')
        ylabel('r [1]')
        title('correlation coefficient (parametric)');
        legend(LEG);
        v=axis; axis([xlim,-.2,1]);
	h = nf;

	for k=1:length(nf)
		subplot(nf(k)); 
		hold on; 
		v = axis;
		if isfield(X,'T') && isfield(X.T,'t0') && isfield(X.T,'t')
			ylim = v([3,4])*[1,.9;0,.1];
			h=patch(X.T.t(X.TC([1,1,end,end,1])),ylim([1,2,2,1,1]),[1,1,1]*0);
%			set(h,'FaceAlpha',.2,'EdgeAlpha',1)
		else	
			plot(X.T.t([X.TC;X.TC])-X.T.t0,[.9;1]*v(4)*ones(1,length(X.TC)),'k-');
		end;
		hold off; 
	end;


elseif strcmp(X.datatype,'TSD_BCI9')  
        if nargin<2,
                clf;
                for k=1:6,
                        nf(k)=subplot(3,2,k);
                end;
        else
                nf=arg2;
        end;

        fid = 1;
        tix = X.tix(1);
        N = length(X.CL);
        c = (N-1)/N;

	fprintf(fid,'Classifier:	%s\n',X.datatype);
	if isfield(X,'hyperparameters')
		f = getfields(X.hyperparameters); 
		for k=1:length(f),
			fprintf(fid,'	%s = %f\n',f{k},getfield(X.hyperparameters,f));
		end; 	
	end; 	
        fprintf(fid,'Error:		     %4.1f %% \n',100-X.ACC00(tix)*100);
        fprintf(fid,'Accuracy:		  %4.1f %% \nspecific Accuracy:	 ',X.ACC00(tix)*100);
        [kap,sd,H,z,OA,SA,MI] = kappa(X.optCMX);
        fprintf(fid,'%4.1f   ',SA*100);
        fprintf(fid,'\nKappa:		     %4.2f %c %4.2f\n',X.KAP00(tix),177,X.Ksd00(tix));
        fprintf(fid,'I(Wolpaw):		 %4.2f bit\n',wolpaw_entropy(X.ACC00(tix),N));
        fprintf(fid,'I(Nykopp):		 %4.2f bit\n',X.I_Nykopp(tix));
        if isfield(X,'I'),
	        fprintf(fid,'I(Continous):	      SUM = %4.2f  [ ',sumskipnan(X.I(tix,:))*c);
        	fprintf(fid,'%4.2f   ',X.I(tix,:));
        	t = X.T; t(t<.5)=NaN;
        	fprintf(fid,' ] \nSTMI:		      %4.2f  [ ',max([sum(X.I,2)]./t)*c);
        	fprintf(fid,'   %4.2f',max(X.I./[t(:,ones(1,size(X.I,2)))-3]));
        	fprintf(fid,' ] \nSNR:		       ')
        	fprintf(fid,'%4.2f   ',X.SNR(tix,:));
        	fprintf(fid,'\ncorrelation (parametric):  ')
        	fprintf(fid,'%4.2f   ',X.r(tix,:));
        	fprintf(fid,'\nrank correlation:	  ');
	        fprintf(fid,'\nAUC:		       ');
	        fprintf(fid,'%4.2f   ',X.AUC(tix,:));
	        fprintf(fid,'\n');
        end; 
        if isfield(X,'rankcorrelation');
                fprintf(fid,'%4.2f   ',X.rankcorrelation(tix,:));
        end;

        xlim = [min(0,X.T(1)),max(X.T(end))];
        
        subplot(nf(1));
        if ~isfield(X,'Labels')
                for k = 1:length(X.CL), Labels{k}=int2str(X.CL(k)); end;
        else
                Labels = X.Labels;
        end;
        %plota(X);
        plot(X.T, 100*[X.ACC00, X.KAP00])
        axis([xlim,-20,100])
        xlabel('time [s]')
        title('Accuracy and Kappa')
        legend({'Accuracy [%]','Kappa [%]'})

        subplot(nf(3));
        plot(X.T,[sum(X.I,2)*c,X.I_Nykopp(:),wolpaw_entropy(X.ACC00,N)])
        legend({'I_{continous}','I_{Nykopp}','I_{Wolpaw}'})
        title('Mutual information')
        ylabel('I [bit]');
        xlabel('time [s]')
        v=axis; axis([xlim,0,1.5]);

        if 0,
                subplot(nf(5));
                plot(X.T,X.MEAN2)
                hold on
                plot(X.T,X.MEAN1)
                hold off
                Title('average output of one-vs-rest classifiers')
                ylabel('TSD')
                xlabel('time [s]')
                tmp = strvcat(Labels);
                legend(cellstr([tmp,repmat('(+)',size(tmp,1),1);tmp,repmat('(-)',size(tmp,1),1)]))
                v=axis; axis([xlim,v(3:4)]);
        else
                subplot(nf(5));
                t = X.T;
                %t(t < 3.5)=NaN;
                t(t < 0.5)=NaN;
                plot(X.T,[sum(X.I,2)*c,X.I_Nykopp(:),wolpaw_entropy(X.ACC00,N)]./repmat(t,1,3))
                legend({'STMI_{C}','STMI_{N}','STMI_{W}'})
                title('Steepness of Mutual information')
                ylabel('STMI [bit/s]');
                xlabel('time [s]')
                v=axis; axis([xlim,0,1]);
        end;

        subplot(nf(2));
        plot(X.T,X.AUC)
        xlabel('time [s]')
        ylabel('AUC [1]')
        title('area-under-the-(ROC)-curve');
        legend(Labels)
        v=axis; axis([xlim,.4,1]);

        subplot(nf(4));
        plot(X.T,[X.I,sum(X.I,2)*c,])
        xlabel('time [s]')
        ylabel('I [bit]');
        title('Mutual Information of continous output');
        legend([Labels,{'SUM'}])
        v=axis; axis([xlim,0,1.5]);

        subplot(nf(6));
        if isfield(X,'N')
        	% show significance interval 
	        alpha = .05;
 	       	ci = tanh(sqrt(2)*erfinv(1-2*alpha)./sqrt(X.N-3));		% confidence interval for alpha of z
	        plot(X.T,X.r,'-',X.T,ci*[-1,1],'k:');
	        LEG = [Labels,{sprintf('alpha=%f',alpha)}];
	else
	        plot(X.T,X.r,'-');
	        LEG = Labels;
        end; 
        xlabel('time [s]')
        ylabel('r [1]')
        title('correlation coefficient (parametric)');
        legend(LEG);
        v=axis; axis([xlim,-.2,1]);
	h = nf;
	
        
elseif strcmp(X.datatype,'TSD_BCI8')    % obsolote
        if ~isfield(X,'T');
                X.T = [1:size(X.ACC00,1)]';
        end;
        h = plot(X.T, 100*[X.ACC00, X.KAP00*[1,1,1,0] + X.Ksd00*[0,-1,1,1]],'-k');
        v = axis; v(3:4)=[-20,129];axis(v)
        set(h(1),'color',[0,0,1]);
        set(h(2),'color',[0,.5,0]);
        set(h(3),'color',[0,1,0]);
        set(h(4),'color',[0,1,0]);
        legend('Accuracy [%]','kappa ± s.d. [%]')
        xlabel('t [s]');

elseif strcmp(X.datatype,'TSD_BCI7')    % obsolete
        if (nargin>1) & strcmpi(arg2,'balken2');
                %elseif strcmpi(X.datatype,'Balken_Diagramm'),

                dy = 0.3;
                F  = 25;
                mn = rs([X.MEAN2(:),X.MEAN1(:)],F,1)';
                pc = rs([X.BCG2(:),X.BCG1(:)],F,1)';
                barh(rs(X.T,F,1),rs([1-2*X.BCG1,X.BCG2*2-1],F,1))

        elseif (nargin>1) & strcmpi(arg2,'balken');
                %elseif strcmpi(X.datatype,'Balken_Diagramm'),

                dy = 0.3;
                F  = 125;
                if isfield(X,'Fs'),
                        if X.Fs==125,
                                F = 50;
                        elseif X.Fs == 128,
                                F = 32;
                        end
                end;
                mn = [rs(X.MEAN2(:),F,1)';rs(X.MEAN1(:),F,1)'];
                pc = [rs(X.BCG2(:),F,1)';rs(X.BCG1(:),F,1)'];
                %barh(rs(X.T,F,1),rs([X.BCG1,X.BCG2],F,1),1)

                samples = (1:length(mn))';
                time = X.T(F:F:end)';
                %time = time - time(1);

                % the following sequence is from R. Scherer
                for k = 1:size(mn, 2),
                        if abs(mn(1,k)) > abs(mn(2,k))
                                patch([0 0 mn(1,k) mn(1,k) 0], [k-dy k+dy k+dy k-dy k-dy], 'k');
                                patch([0 0 mn(2,k) mn(2,k) 0], [k-dy k+dy k+dy k-dy k-dy], 'w');
                        else
                                patch([0 0 mn(2,k) mn(2,k) 0], [k-dy k+dy k+dy k-dy k-dy], 'w');
                                patch([0 0 mn(1,k) mn(1,k) 0], [k-dy k+dy k+dy k-dy k-dy], 'k');
                        end
                        pos1 = min(min(mn(:,k)), 0);
                        pos2 = max(max(mn(:,k)), 0);

                        text(pos1, k, [' ' num2str(100*pc(1,k), '%.0f') '% '], 'HorizontalAlignment', 'Right');
                        text(pos2, k, [' ' num2str(100*pc(2,k), '%.0f') '% '], 'HorizontalAlignment', 'Left');
                end

                mx = max(max(mn));
                mn = min(min(mn));

                %xlim(1.2 * [mn mx]);
                xlim([-1,1]);
                set(gca, 'ytick', 1:length(time), 'yticklabel', num2str(time(:), '%02.2f'));
                ylabel('time in seconds');
                xlabel('distance')
                ylim([0.5 length(samples)+0.5]);
                %	title(caption);

                pos = get(gca, 'position');
                axes('position', pos, 'color', 'none', 'YTick', [], 'XTick', []);
                set(gca, 'YAxisLocation', 'right');

                xlim(1.2 * [mn mx]);
                ylim([0.5 length(samples)+0.5]);
                set(gca, 'ytick', 1:length(time), 'yticklabel', num2str(100*mean(pc)', '%.0f'));

                ylabel('overall classification result')

                %elseif strcmp(X.datatype,'TSD_BCI7') | strcmp(X.datatype,'SNR'), ,
        else
                if nargin<2,
                        clf;
                        for k=1:3,
                                nf(k)=subplot(1,3,k);
                        end;
                else
                        nf=arg2;
                end;
                if isfield(X,'KAP00');
                        if nargin<2,
                                hf(1) = subplot(121);
                                hf(2) = subplot(122);
                        else
                                hf = arg2;
                        end;
                        if ~isfield(X,'T')
                                X.T=(1:size(X.KAP00,1));
                        end;
                        subplot(hf(1));
                        plot(X.T,[X.KAP00,X.ACC00]*100);
                        grid on;
                        v=axis; v(3:4)=[-10,100]; axis(v);
                        ylabel('Kappa [%], Accuracy [%]')
                        xlabel('time t [s]');
                        legend('Kappa', 'Accuracy');

                        subplot(hf(2));
                        h=plot(X.T,[X.I0,X.I]);
                        grid on;
                        v=axis; v(3:4)=[0,1]; axis(v);
                        set(h(end),'linewidQh',2)
                        ylabel('MI [bit]');
                        xlabel('time t [s]');
                        LEG=[];
                        for k = 1:length(X.MD),
                                LEG{k} = int2str(k);
                        end;
                        LEG{k+1}='all';
                        %legend(LEG);
                        return;
                end;

                if ~isfield(X,'Fs'),
                        X.Fs=128;
                end;
                N=length(X.I);
                if isfield(X,'T')
                        t=X.T;%(min(X.T(:)):max(X.T(:)))/X.Fs;
                else
                        t=(1:N)/X.Fs;
                end;

                subplot(nf(1));
                plot(t,X.ERR*100);
                grid on;
                ylabel('Error rate [%]')
                v=axis;v=[0,max(t),0,100];axis(v);

                subplot(nf(2));
                %if strcmp(X.datatype,'SNR'),
                if isfield(X,'MEAN1'),
                        h=plot(t,[X.MEAN1(:),X.MEAN2(:),X.SD1(:),X.SD2(:)]*[1,0,0,0; 0,1,0,0; 1,0,1,0; 1,0,-1,0; 0,1,0,1; 0,1,0,-1]','b',t([1,length(t)]),[0,0],'k');
                else
                        h=plot(t,[X.M1(:),X.M2(:),X.SD1(:),X.SD2(:)]*[1,0,0,0; 0,1,0,0; 1,0,1,0; 1,0,-1,0; 0,1,0,1; 0,1,0,-1]','b',t([1,length(t)]),[0,0],'k');
                end;
                set(h(1),'linewidQh',2);
                set(h(2),'linewidQh',2);
                set(h(7),'linewidQh',2);
                set(h(2),'color','g');
                set(h(5),'color','g');
                set(h(6),'color','g');
                ylabel('Average TSD');
                v=axis;v(1:2)=[0,max(t)];axis(v);

                subplot(nf(3));
                if isfield(X,'T')
                        tmp = repmat(NaN,size(X.T));
                        tmp((X.T>=3.5) & X.ERR<=.50) = 0;
                        plot(t,[X.I,X.I./(X.T-3)+tmp]);
                else
                        plot(t,X.I);
                end;
                ylabel('Mutual Information [bits]')
                v=axis;v=[0,max(t),0,.5];axis(v);

                if length(nf)>3,
                        axes(nf(4))
                        plota(X,'balken');
                elseif 0,
                        plot(t,X.corcov)
                        v=axis;v=[0,10,-1,1];axis(v);
                        grid('on')
                        ylabel('correlation coefficient')
                end;
        end;


elseif strcmp(X.datatype,'ERDS'),
        if nargin==1,
                plotaERDS(X)
        elseif nargin==2,
                plotaERDS(X,arg2)
        end;


elseif strcmp(X.datatype,'QRS_events'),
        ix0 = find(X.EVENT.TYP==hex2dec('0501'));
        if ~isfield(X.EVENT,'CHN')
        	CHAN = 0; 
        elseif ~isempty(ix0)
                CHAN = unique(X.EVENT.CHN(ix0));
        end
        if (length(CHAN)==1)
                ix = X.EVENT.POS(ix0) / X.EVENT.SampleRate;
                if nargin<2,
                        semilogy((ix(1:end-1)+ix(2:end))/2,diff(ix));
                elseif strcmp(arg2,'lin')
                        plot((ix(1:end-1)+ix(2:end))/2,diff(ix));
                else
                        semilogy((ix(1:end-1)+ix(2:end))/2,diff(ix),arg2);
                end;
                ylabel('Inter-Beat-Interval (IBI) [s]')
                xlabel('time [s]')
        else
                for k = 1:length(CHAN),
                        ix = find(X.EVENT.CHN(ix0)==CHAN(k));
                        ix = X.EVENT.POS(ix0(ix)) / X.EVENT.SampleRate;
                        semilogy((ix(1:end-1)+ix(2:end))/2,diff(ix));
                        hold on;
                end;
                hold off;
        end;

        
elseif strcmp(X.datatype,'AMARMA')
        if X.MOP(3)~=0, return; end;
        if (X.MOP(1)==1) & (size(X.AAR,2)==X.MOP(2)+1),
                m0 = X.AAR(:,1)./(1-sum(X.AAR(:,2:end),2));
                ix_aar = 2:X.MOP(2)+1;
        elseif (X.MOP(1)==0) & (size(X.AAR,2)==X.MOP(2)),
                m0 = zeros(size(X.AAR,1),1);
                ix_aar = 1:X.MOP(2);
        else
                return;
        end;

        MODE = [];
        if ~isempty(findstr(upper(arg3),'TIME')),	MODE = [MODE,1]; end;
        if ~isempty(findstr(upper(arg3),'VLHF')),	MODE = [MODE,2]; end;
        if ~isempty(findstr(upper(arg3),'IMAGE')),	MODE = [MODE,3]; end;
        if ~isempty(findstr(upper(arg3),'3D')),		MODE = [MODE,4]; end;
        if ~isempty(findstr(upper(arg3),'n.u.')),	MODE = [MODE,5]; end;
        if ~isempty(findstr(upper(arg3),'LF-BW')),	MODE = [MODE,6]; end;
        if ~isempty(findstr(upper(arg3),'F0')),		MODE = [MODE,7]; end;
        if ~isempty(findstr(upper(arg3),'F0+-B')),	MODE = [MODE,8]; end;
        if ~isempty(findstr(upper(arg3),'ALL')),	MODE = [1,2,3,6]; end;

        if nargin<4,
                if any(MODE==4)
                        hf = gca;
                else
                        for k = 1:length(MODE);
                                hf(k) = subplot(length(MODE),1,k);
                        end;
                end;
        else
                hf = arg4;
        end;
        if ~isfield(X,'T')
                X.T = (0:size(X.AAR,1)-1)';
                X.xlabel = 'beats [1]';
        end;

        HHHH = [];
        K = 0;
        if any(MODE==1)
                K = K + 1;
                subplot(hf(K));
                %plot(X.T,[signal,m0,sqrt([tmp(:,8),X.PE])]);

                if isfield(X,'S')
                        ha=plot(X.T,[X.S,m0,sqrt(X.PE)]);
                        MM1 = max([X.S,m0,sqrt(X.PE)]);
                        MM2 = min([X.S,m0,sqrt(X.PE)]);
                else
                        ha=plot(X.T,[m0,sqrt(X.PE)]);
                        MM1 = max([m0,sqrt(X.PE)]);
                        MM2 = min([m0,sqrt(X.PE)]);
                end;
                MM = [max(MM1),min(MM2)];
                if isfield(X,'EVENT')
                        hold on
                        %text(X.EVENT.POS,repmat(sqrt(max(X.PE))*2,length(X.EVENT.POS),1),'r+');
                        plot(X.T(X.EVENT.POS)*[1,1],MM*[1.2,0;-.20,1],':k');
                        for k = 1:length(X.EVENT.POS),
                                ha(k)=text(X.T(X.EVENT.POS(k)),MM*[1;0],X.EVENT.Desc{k});
                                set(ha(k),'VerticalAlignment','top');
                                set(ha(k),'Rotation',90);
                        end;
                        hold off
                end;
                HHHH = ha;

                v = axis; v(2) = max(X.T); v(3)=0; axis(v);
                ylabel([X.Label,' [',deblank(X.PhysDim),']']);
                hc= colorbar;
                pos=get(gca,'position');
                delete(hc);
                set(gca,'position',pos);
                if isfield(X,'S')
                        legend('Raw','Mean','RMS')
                else
                        legend('mean','RMS')
                end;
        end;

        if prod(size(arg2))==1,
                f0 = repmat(arg2,size(X.AAR,1),1);
        elseif prod(size(arg2))>1,
                f0 = arg2;
        end;

        if any(MODE==2)
                [w,A,B,R,P,F,ip] = ar_spa(X.AAR(:,ix_aar),f0,X.PE);
                ix = (imag(F)==0);

        elseif any(MODE==5) | any(MODE==6) | any(MODE==7) | any(MODE==8)
                [w,A,B,R,P,F,ip] = ar_spa(X.AAR(:,ix_aar),1,X.PE);
                ix = (imag(F)==0);

        end;

        if any(MODE==2)
                ixVLF = ((w>=0)  & (w<.04)); F1 = real(F); F1(~ixVLF)= NaN;
                ixLF  = (w>=.04) & (w<=.15); F2 = real(F); F2(~ixLF) = NaN;
                ixHF  = (w>.15)  & (w<=.4) ; F3 = real(F); F3(~ixHF) = NaN;

                tmp = [sumskipnan(real(F),2), sumskipnan(F1,2), sumskipnan(F2,2), sumskipnan(F3,2)];
                tmp(:,5) = tmp(:,3)./tmp(:,4);
                tmp(:,6) = tmp(:,3)./(tmp(:,1)-tmp(:,2));
                tmp(:,7) = tmp(:,4)./(tmp(:,1)-tmp(:,2));
                tmp(:,8) = sum(tmp(:,2:4),2);

                K = K + 1;
                subplot(hf(K));
                semilogy(X.T,tmp(:,[3,4,8,1])); % 1
                v = axis; v(2) = max(X.T); axis(v);
                ylabel(sprintf('%s [%s^2]',X.Label,X.PhysDim));
                hc= colorbar;
                pos=get(gca,'position');
                delete(hc);
                set(gca,'position',pos);

                %ylabel(sprintf('%s [%s^2/%s]',X.Label,X.PhysDim,'s'));
                legend({'LF','HF','VLF+LF+HF','total'})
        end;

        if any(MODE==5)
                ixVLF = ((w>=0)  & (w<.04)); F1 = real(F); F1(~ixVLF)= NaN;
                ixLF  = (w>=.04) & (w<=.15); F2 = real(F); F2(~ixLF) = NaN;
                %ixHF  = (w>.15)  & (w<=.4) ;
                ixHF  = (w>.17*f0(:,ones(1,size(w,2)))) & (w<=.40*f0(:,ones(1,size(w,2))));
                F3 = real(F); F3(~ixHF) = NaN;

                tmp = [sumskipnan(real(F),2), sumskipnan(F1,2), sumskipnan(F2,2), sumskipnan(F3,2)];
                tmp(:,5) = tmp(:,3)./tmp(:,4);
                tmp(:,6) = tmp(:,3)./(tmp(:,1)-tmp(:,2));
                tmp(:,7) = tmp(:,4)./(tmp(:,1)-tmp(:,2));
                tmp(:,8) = sum(tmp(:,2:4),2);


                K = K + 1;
                subplot(hf(K));
                %semilogy(X.T,tmp(:,[3,4,8])); % 1
                plot(X.T,tmp(:,[6,7])*100); % 1
                v = axis; v(2:4) = [max(X.T),0,100]; axis(v);
                ylabel(sprintf('%s [n.u. %%]',X.Label));
                hc= colorbar;
                pos=get(gca,'position');
                delete(hc);
                set(gca,'position',pos);


                %ylabel(sprintf('%s [%s^2/%s]',X.Label,X.PhysDim,'s'));
                legend({'LF (0.05-0.15Hz)','HF (0.15*f0-0.40*f0)'})
        end;

        if any(MODE==3),
                DN = max(1,ceil(size(X.AAR,1)/1000)); %assume 1000 pixels
                N  = size(X.AAR,1);
                clear h h2 F3
                h2 = repmat(NaN,101,length(1:DN:N));
                t  = repmat(NaN,size(h2,1),1);
                for l = 1:DN:N,  %N/2; [k,size(sdf),N],%length(AR);
                        k = ceil(l/DN);
                        % [h(:,k),F(:,k)] = freqz(sqrt(X.PE(l)/(2*pi*X.AAR(l,1))),[1, -X.AAR(l,2:end)]',128,f0(l,1));
                        %[h2(:,k),F3] = freqz(sqrt(X.PE(l)/(2*pi*X.AAR(l,1))),[1, -X.AAR(l,ix_aar)]',[0:100]'/64,f0(l,1));
                        [h2(:,k),F3] = freqz(sqrt(X.PE(l)/(2*pi*f0(l,1))),[1, -X.AAR(l,ix_aar)]',[0:100]'/64,f0(l,1));
                        h2(find(F3>f0(l,1)/2),k)=NaN;
                        t(k) = X.T(l);
                end;
                K = K + 1;
                subplot(hf(K));
                if 0;%FB{1}=='B';
                        HHHH = imagesc(X.T(1:DN:N),F3,2*log10(abs(h2(:,end:-1:1))));
                elseif 0;
                        HHHH = imagesc(X.T(1:DN:N),F3,2*log10(abs(h2)));
                else
                        HHHH = imagesc(t(~isnan(t)),F3,2*log10(abs(h2(:,~isnan(t)))));
                end;

%                xlabel(X.xlabel);
                ylabel('f [1/s]');
                pos0 = get(gca,'position');

                hc= colorbar;
                %pos1 = get(gca,'position')
                v = get(hc,'yticklabel');
                v = 10.^str2num(v)*2;
                pos2 = get(hc,'position');
                %pos2(1) = 1 - (1 - pos0(1) - pos0(3))/4;
                %pos2(1) = pos0(1) + pos0(3) + pos2(3)/2;
                %pos2(3) = pos2(3)/2
                %set(hc,'yticklabel',v,'position',pos2);

                if 0,isfield(X,'EVENT')
                        hold on
                        %text(X.EVENT.POS,repmat(sqrt(max(X.PE))*2,length(X.EVENT.POS),1),'r+');
                        plot(X.T(X.EVENT.POS)*[1,1],[1.6,0],':k');
                        for k = 1:length(X.EVENT.POS),
                                %[tmp,ix] = min(abs(t-X.EVENT.POS(k)));
                                %plot(ix*[1,1],[1.6,0],':k');
                                %plot(X.T(X.EVENT.POS(k))*[1,1],[1.6,0],':k');
                                ha(k)=text(X.T(X.EVENT.POS(k)),0,X.EVENT.Desc{k});
                                %ha(k)=text(ix,0,X.EVENT.Desc{k});
                                set(ha(k),'VerticalAlignment','top');
                                set(ha(k),'Rotation',90);
                        end;
                        hold off
                end;

                %set(gca,'position',pos0);
                title(sprintf('%s [%s^2/%s]',deblank(X.Label),deblank(X.PhysDim),'s'));
        end;

        if any(MODE==4),
                DN = 8; %10;
                N  = size(X.AAR,1);
                clear h h2 F3
                for l = 1:DN:N,  %N/2; [k,size(sdf),N],%length(AR);
                        k = ceil(l/DN);
                        [h(:,k),F4(:,k)] = freqz(sqrt(X.PE(l)/(2*pi*X.AAR(l,1))),[1, -X.AAR(l,2:end)]',128,f0(l,1));
                        % [h2(:,k),F3] = freqz(sqrt(X.PE(l)/(2*pi*X.AAR(l,1))),[1, -X.AAR(l,2:end)]',[0:70]'/64,f0(l,1));
                        % h2(find(F3>f0(l,1)/2),k)=NaN;
                end;

                K = K + 1;
                ah=subplot(hf(K));
                plot3(repmat(1:DN:ceil(N),128,1),F4,6+2*log10(abs(h)));
                xlabel(X.xlabel);
                ylabel('f [1/s]');
                zlabel(sprintf('%s [%s^2/%s]',X.Label,X.PhysDim,'s'));
                v = get(ah,'zticklabel');
                v = 10.^str2num(v)*2;
                set(ah,'zticklabel',v);

                view([60,-45]);
                pos=get(ah,'position'); %pos(4)=.5;
                set(ah,'position',pos);
                set(1,'PaperUnits','inches','PaperOrientation','portrait','PaperPosition',[0.25 .5 8 10]);
                HHHH = ah;
        end;

        if any(MODE==6)
                ixVLF = ((w>=0)  & (w<.04)); F1 = real(F); F1(~ixVLF)= NaN;
                ixLF  = (w>=.04) & (w<=.15); F2 = real(F); F2(~ixLF) = NaN;
                ixHF  = (w>.15)  & (w<=.4) ; F3 = real(F); F3(~ixHF) = NaN;

                B2 = B; B2(~ixLF) = NaN;
                w2 = w; w2(~ixLF) = NaN;

                HHHH = mean(mean(B2,2).*f0);
                fprintf(1,'Mean LF BandwidQh %f\n',HHHH);

                K = K + 1;
                subplot(hf(K));
                plot(X.T, [mean(w2,2), mean(B2,2).*f0]);
                % hold on
                % errorbar(X.T,mean(w2,2),mean(B2,2));
                % hold off;

                v = axis; v(2:4) = [max(X.T),0,.2]; axis(v);

                hc= colorbar;
                pos=get(gca,'position');
                delete(hc);
                set(gca,'position',pos);

                %legend('f0(LF)','bandwidQh(LF)');
                xlabel(X.xlabel);
                ylabel('f [Hz]');
                legend('f0(LF)','bandwidQh(LF)');
                drawnow

        end;

        if any(MODE==7)
                clim = [0,max(A(:))];
                cm = colormap;

                B2 = B; B2(w<0) = NaN;
                w2 = w; w2(w<0) = NaN;
                for k2 = 1:size(w2,2),
                        for k1 = 1:size(w2,1),
                                if w2(k1,k2)>0,
                                        H = plot(X.T(k1),w2(k1,k2),'d');
                                        set(H,'color',cm(ceil(A(k1,k2)/clim(2)*size(cm,1)),:));
                                end;
                        end;
                end;
        end;
        if any(MODE==8)
                B2 = B; B2(w<0) = NaN;
                w2 = w; w2(w<0) = NaN;
                HHHH = plot(X.T,w2,'d', X.T,w2+B,'v', X.T,w2-B,'^');
        end;
        h = HHHH;


elseif strcmp(X.datatype,'ELPOS2'),
        X = leadidcodexyz(X);
        XYZ = X.ELEC.XYZ; 
        Label = cellstr(X.Label); 
        ix = find(X.LeadIdCode>=0)';
        x = XYZ(:,1); 
        y = XYZ(:,2); 
        R = sqrt(sum(XYZ.^2,2));
        Theta = acos(XYZ(:,3)./R)./sqrt(sum(XYZ(:,1:2).^2,2)); 
        x = x.*Theta;
        y = y.*Theta;
        plot(x,y,'x',x(ix),y(ix),'ro');
        for k=1:length(Label); %ix(:)',
                text(x(k),y(k),[int2str(k),': ',Label{k}]);
        end;
        set(gca,'xtick',0,'ytick',0)
        h = X; 
        if nargin>2,
        	c = arg3;
        	c = ceil(c*10); 
        	for k1=1:size(c,1),
        	for k2=1:size(c,2),
        		%%%% FIXME: Add 2D or 3D connections %%
        		if c(k1,k2)>0,
	        		e= line(x([k1,k2]),y([k1,k2]));
        			set(e,'linewidth',c(k1,k2)); 
        		end;	
        	end;
        	end;
        end;	
               

elseif strncmp(X.datatype,'ELPOS',5),
        X = leadidcodexyz(X);
        XYZ = X.ELEC.XYZ; 
        Label = cellstr(X.Label); 
        ix = find(X.LeadIdCode>0)';
        plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'x',XYZ(ix,1),XYZ(ix,2),XYZ(ix,3),'ro');
        for k=1:length(Label); %ix(:)',
                text(XYZ(k,1),XYZ(k,2),XYZ(k,3),[int2str(k),': ',Label{k}]);
        end;
        set(gca,'xtick',0,'ytick',0)
        h = X; 
        if nargin>2,
        	c = arg3;
        	c = ceil(c*10); 
        	for k1=1:size(c,1),
        	for k2=1:size(c,2),
        		%%%% FIXME: Add 2D or 3D connections %%
        		if c(k1,k2)>0,
	        		e= line(XYZ([k1,k2],1),XYZ([k1,k2],2),XYZ([k1,k2],3));
        			set(e,'linewidth',c(k1,k2)); 
        		end;	
        	end;
        	end;
        end;	
        
        
elseif strcmp(X.datatype,'REV'),
        if nargin<2
                arg2=[];
        end;

        if strcmpi(arg2,'3d')
                R=X.REV;
                [mxR,ix]=min(R);
                [myR,iy]=min(R');
                [IX,IY]=find(R==min(R(:)));
                [IX,IY],
                h=mesh(R');hold on
                plot3(ix,1:length(ix),mxR,'go-')
                plot3(1:length(iy),iy,myR,'cx-')
                plot3(IX,IY,min(R(:)),'r*');
                hold off
                %v=axis;v(5:6)=[0.1,.2];axis(v);
                a=gca;
                set(a,'xtick',X.P);
                t=1:5:length(X.UC);
                set(a,'ytick',t);
                set(a,'yticklabel',X.UC(t));
                ylabel('update coefficent UC')
                xlabel('model order p')
                zlabel('REV = MSE/MSY')
        else  %if strcmpi(arg2,'2d')
                subplot(121)
                semilogx(X.UC,X.REV')
                v=axis;v=[1e-6,1e-1,0.1,.15];axis(v);
                xlabel('update coefficent UC')
                ylabel('REV = MSE/MSY')

                subplot(122)
                h=plot(X.P,X.REV)
                v=axis;v(3:4)=[0.1,.15];axis(v);
                xlabel('model order p')
                ylabel('REV = MSE/MSY')

        end;

elseif strcmp(X.datatype,'MDist-matrix'),
        if nargin<2
                arg2=[];
        end;
        h=mesh(X.D)
        h=get(gca,'xtick');
        %set(gca,'xticklabel',X.t(h(h>0)))
        h=get(gca,'ytick');
        %set(gca,'yticklabel',X.t(h(h>0)))

elseif strcmp(X.datatype,'MD'),
        if nargin<2
                arg2=[];
        end;

        if strcmpi(arg2,'3d')
                R=X.REV;
                [mxR,ix]=max(R);
                [myR,iy]=max(R');
                [IX,IY]=find(R==max(R(:)));
                h=mesh(R');hold on
                plot3(ix,1:length(ix),mxR,'go-')
                plot3(1:length(iy),iy,myR,'cx-')
                plot3(IX,IY,min(R(:)),'r*');
                hold off
                %		v=axis;v(5:6)=[0.1,.2];axis(v);
                a=gca;
                set(a,'xtick',X.P);
                t=1:5:length(X.UC);
                set(a,'ytick',t);
                set(a,'yticklabel',X.UC(t));
                ylabel('update coefficent UC')
                xlabel('model order p')
                zlabel('log_{10}(MD)')
        else  %if strcmpi(arg2,'2d')
                subplot(121)
                semilogx(X.UC,X.REV')
                %		v=axis;v=[1e-6,1e-1,0.1,.15];axis(v);
                xlabel('update coefficent UC')
                ylabel('log_{10}(MD)')

                subplot(122)
                h=plot(X.P,X.REV)
                %	       v=axis;v(3:4)=[0.1,.15];axis(v);
                xlabel('model order p')
                ylabel('log_{10}(MD)')
        end;
end;

if nargout,
        H = h;
end;
