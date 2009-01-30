function [OCO]=oco(fn,eegchan,trigchan,Mode)
% OCO - evaluates online classification output 
%
% CC = oco(fn)
%  
% CC = oco(...,Mode)
%       Mode 'E' provides extended analysis
%  
%  e.g. 
%  CC=OCO('x21fb*')
%  CC=OCO({'x21fb1.mat','x21fb2.mat'})
%  CC=OCO(fn,[1,3],4,128)
%
% default: 
% 	eegchan=[1,3];
% 	trigchan=4;
% 	Fs=128;


%	$Revision: 1.1 $
%	$Id: oco.m,v 1.1 2009-01-30 06:04:50 arno Exp $
%	Copyright (C) 1997-2004 by Alois Schloegl 
%	a.schloegl@ieee.org	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 2 of the License, or (at your option) any later version.
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


cname=computer;
if (nargin<2),
        eegchan = nan;
elseif ischar(eegchan)
        Mode = eegchan;
        eegchan = nan;
end;
if nargin<3,
        trigchan = 4;
elseif ischar(trigchan)
        Mode = trigchan;
        trigchan = 4;
end;
if exist('Mode')~=1,
        Mode = '';
end

M0  = 7;
MOP = 3;
uc  = 30:5:80;

%lsg1234;
[S,H] = sload(fn,'TSD');
if size(S,2)<=H.NS,
        fprintf(2,'OCO: %s.tsd not available\n',H.FILE(1).Name);
        return;
end;
Fs    = H.SampleRate;

if strcmp(H.TYPE,'GDF'),
        ix0 = (H.EVENT.TYP>hex2dec('300')) & (H.EVENT.TYP<hex2dec('30d')); 
        H.Classlabel = mod(H.EVENT.TYP(ix0),256);

        ix1 = (H.EVENT.TYP==hex2dec('300')); 
        if ~isempty(ix1)
                TRIG = H.EVENT.POS(ix1);
        else
                TRIG = H.EVENT.POS(ix0);
        end
else
        TRIG = gettrigger(S(:,trigchan));
end
cl = H.Classlabel(:);

if length(TRIG)~=length(cl);
        fprintf(2,'number of Triggers (%i) does not fit size of class information (%i)',length(TRIG),length(cl));
        return;        
end;

if isfield(H,'TriggerOffset');
        TRIG = TRIG - H.TriggerOffset*Fs/1000;
        LEN  = (H.BCI.Paradigm.TrialDuration/1000+1)*Fs;
else
        TRIG = TRIG - 2*Fs;
        LEN  = 9*Fs;
end;
OCO = bci4eval(S(:,H.NS+1),TRIG, cl, 1, 9*Fs, Fs);

if ~strcmp(Mode,'E'),
%        plota(OCO,'balken');
else
        % Muscle Detection 
        [INI,s,E,arti] = detectmuscle2(S(:,1:3));
        subplot(211);
        if size(S,2)>H.NS,
                plot((1:size(S,1))'/H.SampleRate,[S(:,1:3),E-50,arti-90,S(:,H.NS+1)*10-150]);
        else
                plot((1:size(S,1))'/H.SampleRate,[S,E-50,arti-90]);
        end;
        title(H.FileName)
        for k=5:8, hf(k) = subplot(2,4,k); end;
        plota(OCO,hf(5:8));
end;
