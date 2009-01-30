function HDR=bdf2biosig_events(EVENT)
% BDF2BIOSIG_EVENTS converts BDF Status channel into BioSig Event codes. 
%
%  HDR = bdf2biosig_events(arg1)   
%
%  arg1 can be an HDR-struct containg  HDR.EVENT.Desc 
%  or a struct containing EVENT.Desc 
%  or a cell-array Desc
%  or a char-array Desc 
% 
%  Warning: Approximately 32 (out of 510) events are currently not supported. 
%  For your data, you can check this with this command: 
%  HDR.EVENT.Desc(HDR.EVENT.TYP==0) 
% 
% see also: doc/eventcodes.txt

%	$Id: bdf2biosig_events.m,v 1.1 2009-01-30 06:04:40 arno Exp $
%	Copyright (C) 2007 by Alois Schloegl <a.schloegl@ieee.org>	
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


HDR.EVENT.Desc = []; 
if isfield(EVENT,'EVENT')
	HDR = EVENT; 
elseif isfield(EVENT,'Desc')
	HDR.EVENT = EVENT; 
elseif iscell(EVENT)
	HDR.EVENT.Desc = EVENT; 
elseif ischar(EVENT)
	for k=1:size(EVENT,1),
		HDR.EVENT.Desc{k} = EVENT(k,:); 
	end; 
else 
	error('unknown input argument')
end; 	

status_channel = strmatch('Status',HDR.Label); 

if strcmp(HDR.TYPE,'BDF')

elseif isfield(HDR,'data')
	tmp = strmatch('Status',HDR.Label);
%	HDR.BDF.Status = tmp;
end; 

t = HDR.BDF.ANNONS;

ix1 = diff(double([0;bitand(t,2^16)]));	% start of epoch
ix2 = diff(double([0;bitand(t,2^16-1)]));	% labels 

% defines mapping of the BDF-status channel to BioSig event codes 			
switch 1,    % determines default decoding
case 1,	
	% epoching information is derived from bit17
	% only lower 8 bits are supported
	POS = [find(ix1>0);find(ix2>0);find(ix1<0);find(ix2<0)];
	TYP = [repmat(hex2dec('7ffe'),sum(ix1>0),1); bitand(t(ix2>0),255); repmat(hex2dec('fffe'),sum(ix1<0),1); bitor(bitand(t(find(ix2<0)-1),255),2^15)];

case 2, 
	% suggested decoding if standardized event codes (according to 
	% .../biosig/doc/eventcodes.txt) are used  
	POS = [find(ix2>0)];
	TYP = [bitand(t(ix2>0),2^16-1)];
	
case 3,
	% the status word is decoded only if bit17=0. 
	% raising and falling edges are considered
	t = bitand(HDR.BDF.ANNONS,2^16-1);
	t(~~bitand(HDR.BDF.ANNONS,2^16)) = 0; 
	ix2 = diff(double([0;bitand(t,2^16-1)]));	% labels 
	POS = [find(ix2)];
 	TYP = [t(ix2>0); t(find(ix2<0)-1)+hex2dec('8000')];
		
case 4,
	% the status word is decoded only if bit17=0. 
	% only raising edges are considered
	t = bitand(HDR.BDF.ANNONS,2^16-1);
	t(~~bitand(HDR.BDF.ANNONS,2^16)) = 0; 
	ix2 = diff(double([0;bitand(t,2^16-1)]));	% labels 
	POS = [find(ix2>0)];
	TYP = [t(ix2>0)];
		
case 5,
	% the status word is decoded only if bit17=0. 
	% raising and falling edges are considered
	t = bitand(HDR.BDF.ANNONS,255);			% only bit1-8 are considered, useful if bit9-16 are open/undefined
	t(~~bitand(HDR.BDF.ANNONS,2^16)) = 0; 
	ix2 = diff(double([0;bitand(t,2^16-1)]));	% labels 
	POS = [find(ix2)];
 	TYP = [t(ix2>0); t(find(ix2<0)-1)+hex2dec('8000')];
		
case 6,
	% the status word is decoded only if bit17=0. 
	% only raising edges are considered
	t = bitand(HDR.BDF.ANNONS,255);			% only bit1-8 are considered, useful if bit9-16 are open/undefined
	t(~~bitand(HDR.BDF.ANNONS,2^16)) = 0; 
	ix2 = diff(double([0;bitand(t,2^16-1)]));	% labels 
	POS = [find(ix2>0)];
	TYP = [t(ix2>0)];
		
case 99,
	% not recommended, because it could break some functionality in BioSig 
	POS = [find(ix2>0);find(ix2<0)];
 	TYP = [bitand(t(ix2>0),2^16-1); bitor(bitand(t(find(ix2<0)-1),2^16-1),2^15)];

end;
[HDR.EVENT.POS,ix] = sort(POS);
HDR.EVENT.TYP = TYP(ix);


%%%% BDF Trigger and status 
t = bitand(HDR.BDF.ANNONS,hex2dec('00ffff')); 
ix = diff(double([0;t]));
HDR.BDF.Trigger.POS = find(ix);
HDR.BDF.Trigger.TYP = t(HDR.BDF.Trigger.POS); 			

t = bitand(bitshift(HDR.BDF.ANNONS,-16),hex2dec('00ff')); 
ix = diff(double([0;t]));
HDR.BDF.Status.POS = find(ix);
HDR.BDF.Status.TYP = t(HDR.BDF.Status.POS); 			

%HDR.BDF.ANNONS = []; 	% not needed anymore, saves memory
