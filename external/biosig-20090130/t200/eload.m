function [HDR] = eload(filename,Fs)
% ELOAD loads EVENT data 
% Event information is often stored in different formats. 
% ELOAD tries to load different formats into a unified 
% form 
% 
% HDR = eload(filename)
%
% filename	Filename of Event information 
% HDR.EVENT contains the EVENT information
% 
% 
% see also: SLOAD, SVIEW, SOPEN 
%


%	$Revision: 1.1 $
%	$Id: eload.m,v 1.1 2009-01-30 06:04:40 arno Exp $
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


HDR = getfiletype(filename);

if strcmp(HDR.TYPE,'MAT')
        tmp = load('-mat',filename);
        if isfield(tmp,'eventmatrix') & isfield(tmp,'samplerate') 
                %%% F. Einspieler's Event information 
                HDR.EVENT.POS = tmp.eventmatrix(:,1);
                HDR.EVENT.TYP = tmp.eventmatrix(:,2);
                HDR.EVENT.CHN = tmp.eventmatrix(:,3);
                HDR.EVENT.DUR = tmp.eventmatrix(:,4);
                HDR.EVENT.Fs  = tmp.samplerate;
                HDR.TYPE = 'EVENT';
                
        elseif isfield(tmp,'EVENT') 
                HDR.EVENT = EVENT; 
                HDR.TYPE = 'EVENT';
        end;
        
elseif strcmp(HDR.TYPE,'GDF');
        H = sopen(HDR,'r'); H=sclose(H);
        HDR.EVENT = H.EVENT; 
        HDR.EVENT.Fs = H.SampleRate; 
        HDR.TYPE = 'EVENT';
        
elseif strncmp(HDR.TYPE,'BrainVision',11);
        HDR = sopen(HDR,'r'); HDR=sclose(HDR); 
	if isfield(HDR.EVENT,'TeegType')
		ix = strmatch('New Segment',HDR.EVENT.TeegType); 
		HDR.EVENT.TYP(ix)=hex2dec('7ffe'); 
	end; 
	for k1 = 1:length(HDR.EVENT.Desc)
		tmp = HDR.EVENT.Desc{k1};
		%HDR.TRIG = HDR.EVENT.POS(HDR.EVENT.TYP<10); 
		if 0,
			
	        elseif strncmp(tmp,'TargetCode',10)
	        	HDR.EVENT.TYP(k1) = str2double(tmp(11:12))+hex2dec('0300'); 
	        elseif strcmp(tmp,'BeginOfTrial')
	        	HDR.EVENT.TYP(k1) = hex2dec('0300'); 
                elseif strcmp(tmp,'hit')
	        	HDR.EVENT.TYP(k1) = hex2dec('0381'); 
	        elseif strcmp(tmp,'wrong')
	        	HDR.EVENT.TYP(k1) = hex2dec('0382'); 

	% eye movements
	        elseif strcmpi(tmp,'augen links')	
	        	HDR.EVENT.TYP(k1) = hex2dec('0431');
	        elseif strcmpi(tmp,'augen rechts')	
	        	HDR.EVENT.TYP(k1) = hex2dec('0432');
	        elseif strcmpi(tmp,'augen hoch') | strcmpi(tmp,'augen oben')		
	        	HDR.EVENT.TYP(k1) = hex2dec('0433');
	        elseif strcmpi(tmp,'augen unten') | strcmpi(tmp,'augen runter')		
	        	HDR.EVENT.TYP(k1) = hex2dec('0434');
	        elseif strcmpi(tmp,'augen offen')	
	        	HDR.EVENT.TYP(k1) = hex2dec('8430');
	        elseif strcmpi(tmp,'augen zu')	
	        	HDR.EVENT.TYP(k1) = hex2dec('0430');
	        elseif strcmp(tmp,'blinzeln')	
	        	HDR.EVENT.TYP(k1) = hex2dec('0439'); 
	
	% muscle movements 
	        elseif strcmp(tmp,'EMG links')
	        	HDR.EVENT.TYP(k1) = hex2dec('0441'); 
	        elseif strcmp(tmp,'EMG rechts')
	        	HDR.EVENT.TYP(k1) = hex2dec('0442'); 
	        elseif strcmpi(tmp,'kopf bewegen')
	        	HDR.EVENT.TYP(k1) = hex2dec('0443'); 
	        elseif strcmp(tmp,'zunge an')
	        	HDR.EVENT.TYP(k1) = hex2dec('0444'); 
	        elseif strcmp(tmp,'Kiefer anspannen')
	        	HDR.EVENT.TYP(k1) = hex2dec('0445'); 
	        elseif strcmp(tmp,'zunge aus')
	        	HDR.EVENT.TYP(k1) = hex2dec('8444'); 
	        elseif strcmp(tmp,'kopf beißen') | strcmp(tmp,'kopf beißen'),
	        	HDR.EVENT.TYP(k1) = hex2dec('0446'); 
	        elseif strcmp(tmp,'EMG fuss')
	        	HDR.EVENT.TYP(k1) = hex2dec('0447'); 
	        elseif strcmp(tmp,'Arme bewegen')
	        	HDR.EVENT.TYP(k1) = hex2dec('0449'); 

	        elseif strncmp(tmp,'S',1)
	        	n = str2double(tmp(2:end)); 
			if n==11,	% hit (left)
			       	HDR.EVENT.TYP(k1) = hex2dec('0381'); 
			elseif n==12,	% hit (right)
		        	HDR.EVENT.TYP(k1) = hex2dec('0381'); 
	        	elseif n==21,	% miss (left)
		        	HDR.EVENT.TYP(k1) = hex2dec('0382'); 
			elseif n==22,	% miss (right)
		        	HDR.EVENT.TYP(k1) = hex2dec('0382'); 
			elseif n==60,	% feedback onset
		        	HDR.EVENT.TYP(k1) = hex2dec('030d'); 
			else
		        	HDR.EVENT.TYP(k1) = n; 
			end; 
        	
	        elseif strcmp(tmp,'s') | strcmp(tmp,'stop') | strcmp(tmp,'stopp'),
	        	HDR.EVENT.TYP(k1) = bitxor(hex2dec('8300'),HDR.EVENT.TYP(k1-1)); 
	        	
	        elseif ~isempty(tmp)
	        	[n,v,s] = str2double(tmp(2:end)); 
	        	if (length(n)==1) & (~v)
	        		HDR.EVENT.TYP(k1) = n; 
	       		end; 
	        end; 	
	end; 
	HDR.EVENT.TYP = HDR.EVENT.TYP(:); 

	if isfield(HDR.EVENT,'POS'); 
	       	ix1 = find(HDR.EVENT.TYP<10); 
	       	ix2 = find(HDR.EVENT.TYP==100); 
		HDR.EVENT.TYP(ix2,1) = HDR.EVENT.TYP(ix2-1)+hex2dec('8000'); 
		ix0 = find((HDR.EVENT.TYP>0)&(HDR.EVENT.TYP<10));
		HDR.TRIG = HDR.EVENT.POS(ix0); 
		HDR.Classlabel = HDR.EVENT.TYP(ix0); 
	end; 

        HDR = bv2biosig_events(H); 
        
        %%% Artifact database of the sleep EEG 
elseif strcmp(HDR.FILE.Ext,'txt') & strmatch(HDR.FILE.Name,['h000201';'h000901';'h001001']);
        HDR.EVENT = adb2event(filename,100);        
        HDR.TYPE = 'EVENT';
elseif strcmp(HDR.FILE.Ext,'txt') & strmatch(HDR.FILE.Name,['b000101';'b000401';'c000701';'c001701';'m000401';'m000901']);
        HDR.EVENT = adb2event(filename,200);        
        HDR.TYPE = 'EVENT';
elseif strcmp(HDR.FILE.Ext,'txt') & strmatch(HDR.FILE.Name,['n000101';'n000401';'p000101';'p000201';'s000201']);
        HDR.EVENT = adb2event(filename,256);        
        HDR.TYPE = 'EVENT';
elseif strcmp(HDR.FILE.Ext,'txt') & strmatch(HDR.FILE.Name,['u000601']);
        HDR.EVENT = adb2event(filename,400);        
        HDR.TYPE = 'EVENT';

elseif strcmp(HDR.TYPE,'WSCORE_EVENT')
        %HDR.EVENT.POS = HDR.EVENT.POS;         % already defined
        HDR.EVENT.TYP = HDR.EVENT.WSCORETYP;         % code assignment not
        fprintf(2,'Warning ELOAD: Event Codes in file %s do not not follow the standard codes of BIOSIG.\n',filename);
        %defined 
        
else
        fprintf(2,'Warning ELOAD: file %s is not recognized as event file.\n',filename);
        
end;