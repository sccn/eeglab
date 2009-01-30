function [argout,s]=sview(s,varargin),
% SVIEW - a simple signal viewer 
%    SVIEW(filename)
%    SVIEW(HDR)
%    SVIEW(S,HDR)
%
% See also: SLOAD 

%	$Id: sview.m,v 1.1 2009-01-30 06:04:51 arno Exp $ 
%	Copyright (c) 2004,2006,2008 by Alois Schlögl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

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

if length(varargin),
        H=varargin{1}; 
        arg2=varargin{1}; 
else
	arg2='';     
	H = [];    
end;

	if ~isempty(H)
		CHAN = 1:size(s,2);

	elseif ischar(s) || iscell(s),
		[s,H] = sload(s,0,'OVERFLOWDETECTION:OFF',varargin{:});
		CHAN = 1:size(s,2);
	
	elseif isstruct(s)
	        [s,H] = sload(s);
		CHAN = 1:size(s,2);
        
	elseif isnumeric(s) & (length(H.InChanSelect)==size(s,2))
	        CHAN = 1:size(s,2);
	        H.Label = H.Label(H.InChanSelect,:);
	else    
        	return;
        end;	

%if strmatch(H.TYPE,{'BMP','PBMA','PGMA','PPMA','PBMB','PGMB','PPMB','XPM'}),
if isfield(H,'IMAGE');
        if exist('OCTAVE_VERSION','builtin')
                if (length(size(s))==3) & (size(s,3)==3)
                        imshow(s(:,:,1),s(:,:,2),s(:,:,3));
		else
		        imshow(s);
                end;	
        else	
                image(s);
        end;
        argout=H;
        return;

elseif ischar(arg2) & (strcmp(H.TYPE,'ELPOS') | (isfield(H,'ELEC') && strncmpi(arg2,'ELPOS',5)));
	XYZ = H.ELEC.XYZ;
	K = 1:size(XYZ,1); 
		
	if ~isempty(strfind(arg2,'2d'))
		T=0:.001:2*pi;

		r = sqrt(sum(XYZ(:,1:2).^2,2));
		Rmax = max(r);
		R = Rmax * 11/10; 
		r = 10/11*atan(r./XYZ(:,3));
		XY = XYZ(:,1:2).*r(:,[1,1]);
		plot(XY(:,1),XY(:,2),'x',XY(K,1),XY(K,2),'ro',R*sin(T),R*cos(T),'b-',-R+R/10*sin(-T/2),R/11*cos(-T/2),'b-',R/11*sin(T/2)+R,R/11*cos(T/2),'b-',[-R 0 R]/11,[R R*1.1 R],'b-');
		for k=1:size(XYZ,1),
			text(XY(k,1),XY(k,2),H.Label(k,:));
		end;
	else
		plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'x',XYZ(K,1),XYZ(K,2),XYZ(K,3),'ro');
		for k=1:size(XYZ,1),
			text(XYZ(k,1),XYZ(k,2),XYZ(k,3),H.Label(k,:));
		end;
	end;
        argout=H;
        return;

        
elseif all((H.LeadIdCode>0) & (H.LeadIdCode<256)) && (H.SPR*H.NRec==H.SampleRate*10) && all(H.PhysDimCode==H.PhysDimCode(1)) 
	% 12-lead, 10s ECG 
	xlen = 5; 
	d = repmat(NaN,[H.SampleRate*xlen,12]);
	pos = [1,2,61:64,3:8];
	x = leadidcodexyz(pos); 
	[tmp, scale] = physicalunits(H.PhysDimCode);
	CH = zeros(1,12);
	for k=1:length(pos),
		ch = find(pos(k)==H.LeadIdCode);
		if ~isempty(ch)
			CH(ch) = k;
			d(1:xlen*H.SampleRate-1,k)=s(1:xlen*H.SampleRate-1,ch);
		else 	
			switch k
			case 3
				d(:,3)=d(:,1:2)*[-1;1];
			case 4
				d(:,4)=d(:,1:2)*[1;1]/2;
			case 5
				d(:,5)=d(:,1:2)*[1;-1/2];
			case 6
				d(:,6)=d(:,1:2)*[-1/2;1];
			end;
		end; 	
	end;

	d = reshape(permute(reshape(d,[H.SampleRate*xlen,3,4]),[1,3,2]),[H.SampleRate*xlen*4,3]);
	d = [(abs([-(H.SampleRate/2):(H.SampleRate/2)])' < (H.SampleRate/4))*ones(1,3)*1e-3/scale(1); repmat(NaN,1,3); d];	
	
	plot([1:H.SampleRate*(xlen*4+1)+2]'/H.SampleRate, d + ones(size(d,1),1)*[0,-1,-2]*3*1e-3/scale(1),'k');
	set(gca,'xMinorTick','on','yMinorTick','on','xMinorGrid','on','yMinorGrid','on');
	grid on
	for k=1:length(pos),
		h=text(xlen*floor((k-1)/3)+1,(1-3*mod(k-1,3))*1e-3/scale(1),x.Label(k));
	end;
	text(.1,-2e-4/scale(1),'0.5s / 1 mV'); 
	set(gca,'xlim',[0,4*xlen+1]);

	gender='UMFX'; sex=gender(H.Patient.Sex+1);
	if isfield(H.Patient,'Age') age=H.Patient.Age;
	elseif isfield(H.Patient,'Birthday') && (datenum(H.Patient.Birthday)>0) 
		age = datenum(H.T0-H.Patient.Birthday)/365.25;
	else age = -1;
	end;		
	id = H.Patient.Id; id(id=='_')=' ';
	title(sprintf('Patient Id: %s   Age=%i y  Sex=%c',id,age,sex));

	tmp = physicalunits(H.PhysDimCode(1));
	ylabel(tmp{1});
	xlabel('time [s]');
	return;

        
elseif strcmp(H.TYPE,'unknown');
        argout=H;
        return;
end;

%s(abs(s)>1e3)=NaN;

[p,f,e]=fileparts(H.FileName);
%fn=dir(fullfile(p,[f,'EOG',e]));
if 1,   % no EOG corrections
elseif 1
                nx = sparse(28,3);
                nx(23:24,1) = [1;-1];
                nx(24:25,2) = [1;-1];
                nx(26,3) = 1;
                [R,s] = regress_eog(H.FileName,1:22,nx);
		size(s),
		[u,s0,v] = svd(s(~any(isnan(s),2),:),0);
		size(u)
		s = s*v;
elseif 1
        [R,s] = regress_eog(s,1:60,61);
end;

if ~isfield(H,'Label'),
        LEG = '';
elseif all(size(CHAN)>1)
        LEG = 'rereferenced';
elseif all(CHAN>0)
        LEG = H.Label(CHAN);
elseif size(H.Label,1)<H.NS,
      	LEG = H.Label;
else
        LEG = H.Label(CHAN);
end;

t = detrend(s); t = t(:); 
t(isnan(t))=median(t);
dd = max(t)-min(t);
dd = 1000; s=center(s); 
dd = max(std(s))*5;
s = zscore(s); dd = 20; % 
H.AS.TIMECHAN = [strmatch('Time',H.Label);strmatch('Zeit',H.Label)];
FLAG.tmp = (length(H.FILE)==1) & (prod(size(H.AS.TIMECHAN))==1);
if FLAG.tmp,
        % this construct is necessary for compatibility with Octave and Matlab 5.3
        FLAG.tmp = any(H.AS.TIMECHAN==H.InChanSelect);
end;

if FLAG.tmp,
        t = s(:,H.AS.TIMECHAN);
        s(:,H.AS.TIMECHAN) = NaN;
        X_Label = [H.Label(H.AS.TIMECHAN,:),' [',H.PhysDim{H.AS.TIMECHAN},']'];
elseif isnan(H.SampleRate)
        t = (1:size(s,1))';
        X_Label = 'samples';
else
        t = (1:size(s,1))'/H.SampleRate;
        X_Label = 'time [s]';
end;

if isfield(H,'PLX'),
elseif isfield(H,'SegLen'), 
        EVENT.POS = H.SegLen(1:end-1)'+1;
        EVENT.DUR = zeros(size(EVENT.POS));
        EVENT.Desc = cellstr(repmat(' ',length(H.FILE),1));
%        t = (1:size(s,1))';
%        X_Label = 'samples';
elseif isfield(H.EVENT,'Desc'),
        EVENT = H.EVENT;
elseif isfield(H.EVENT,'CodeDesc'),
        EVENT = H.EVENT;
        EVENT.CodeDesc{32766} = 'New Segment';
        EVENT.CodeDesc{32767} = 'sparse sample';
	EVENT.Desc = EVENT.CodeDesc(H.EVENT.TYP);	
elseif ~isfield(H.EVENT,'Desc') & (length(H.EVENT.TYP)>0),
	[p,f,e]=fileparts(which('sopen.m'));
	[p,f,e]=fileparts(p);
        g = sopen(fullfile(p,'doc/eventcodes.txt'));g=sclose(g);
        ix = sparse(g.EVENT.CodeIndex,1,1:length(g.EVENT.CodeIndex));
        EVENT = H.EVENT;
        try,
                EVENT.Desc = {g.EVENT.CodeDesc{ix(H.EVENT.TYP)}};
        catch
                fprintf(2,'SVIEW: unknown eventcodes in file %s',H.FileName)
        end
else
        EVENT.POS = [];
end;


if H.NS==1,
        plot(t(:),s,'-');
else
        plot(t(:),((s+(ones(size(s,1),1)*(1:size(s,2)))*dd/(-2)+4*dd)),'-');
end;
%box('off'); 


ix = find(EVENT.POS>0 & EVENT.POS<=length(t));
if isfield(EVENT,'Desc') & ~isempty(ix) %& ~exist('OCTAVE_VERSION','builtin')
	ix2 = find((EVENT.POS>0) & (EVENT.POS<=length(t))); % & (EVENT.DUR>0));
        v = axis;
        hold on;
        N  = length(EVENT.POS);
        plot([1;1]*t(EVENT.POS(ix))',v(3:4),':k');
        if length(EVENT.DUR)==length(EVENT.POS),
        for k=1:length(ix2),
		x = t(EVENT.POS(ix2(k))+[0,EVENT.DUR(ix2(k))])';
		y = v(3:4);
		if exist('OCTAVE_VERSION','builtin');
			patch(x([1,1,2,2,1]),y([1,2,2,1,1]),'b');
		else
			ha=patch(x([1,1,2,2,1]),y([1,2,2,1,1]),[1,1,1]/4*3);
			set(ha,'FaceAlpha',.5);
			set(ha,'EdgeAlpha',.5);
		end;
        end;
        end;
        ix = ix2;
	for k=1:length(ix),
                txt = EVENT.Desc{ix(k)};
                if isempty(txt),txt='';end; 
                txt(txt=='_')=' ';
		if exist('text','builtin');
			if exist('OCTAVE_VERSION','builtin');
				text(t(EVENT.POS(ix(k))),v(3:4)*[-.2;.8],txt);
			else 
			        ha=text(t(EVENT.POS(ix(k))),v(3:4)*[-.2;.8],txt);
        	    		set(ha,'Rotation',90,'VerticalAlignment','top');
        	        end
                end;
        end;	
        %v(4) = v(3:4)*[-.7;1.3]; axis(v);
        hold off;
end;

if 0, length(H.EVENT.POS) > 0,
        hold on;
        if 0, 
        elseif isfield(H.EVENT,'DUR') & isfield(H.EVENT,'CHN');
                plot([H.EVENT.POS,H.EVENT.POS+H.EVENT.DUR]'/H.SampleRate,[dd;dd]*H.EVENT.CHN','+-');    
                
        elseif isfield(H.EVENT,'CHN');
                plot(H.EVENT.POS/H.SampleRate, H.EVENT.CHN*dd, 'x');
                
        elseif isfield(H.EVENT,'DUR');
                plot([H.EVENT.POS,H.EVENT.POS+H.EVENT.DUR]'/H.SampleRate,[dd;dd]*H.EVENT.CHN','+-');    
                
        else; 
                plot(H.EVENT.POS/H.SampleRate,dd*ones(H.EVENT.N,1),'^');    
                
        end;
        hold off;
        LEG = strvcat(LEG,'Events');
end;

if length(H.FILE)==1,
        tmp = H.FileName; tmp(find(tmp==92))='/'; tmp(find(tmp=='_'))=' ';
else
        tmp = '';
end;
title([tmp, ' generated with BIOSIG tools for Octave and Matlab(R)']);
xlabel(X_Label);
PhysDim = '';
if ~isempty(H.PhysDim),
        PhysDim = deblank(H.PhysDim{CHAN(1)});
end;
ylabel(sprintf('Amplitude [%s]',PhysDim));

if exist('OCTAVE_VERSION')<5;
        if ~isempty(LEG);
                legend(LEG);
        end;
end;

if nargout, 
        argout=H;
end;
