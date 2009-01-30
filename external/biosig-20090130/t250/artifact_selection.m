function [HDR,A] = artifact_selection(fn,t1,t2)
% ARTIFACT_SELECTION returns the selected triggered trials with 
% artifacts. The length of the trial is defined by t1 and t2 in 
% seconds.  
%
% [HDR] = artifact_selection(filename,[t1,t2])
% [HDR] = artifact_selection(HDR,[t1,t2])
%    uses EVENT information in filename or HDR.EVENT. 
%
% The artifact selection is available in HDR.ArtifactSelection. 
% HDR.ArtifactSelection is vector with the same length than the list
% of trigger points (event 0x0300). A value of 1 indicates an
% artifact, a value of 0 means the trial is free from any artifacts.
% 	
% In case the trigger information and the artifact scoring is 
% stored in separate files, the information of several files can be
% merged.
%
% [HDR] = artifact_selection({sourcefile,eventfile1,eventfile2,...},[t1,t2])
%
% All files can be defined by their filename, or by the BIOSIG HDR-struct.
% The header of the first file is merged with the Event information of 
%    all other files. 
%
% see also: TLOAD, SLOAD, SOPEN,  

%	$Revision: 1.1 $
% 	$Id: artifact_selection.m,v 1.1 2009-01-30 06:04:42 arno Exp $
%	Copyright (c) 2004-2005,2007 by Alois Schloegl <a.schloegl@ieee.org>
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


% load and merge all artifact information
if ischar(fn)
	HDR = sopen(fn,'r'); HDR = sclose(HDR);
	EVENTMATRIX = [HDR.EVENT.POS,HDR.EVENT.TYP,HDR.EVENT.CHN,HDR.EVENT.DUR];
elseif iscell(fn)
	for k = 1:length(fn)
		if ischar(fn{k})
			h = sopen(fn{k},'r'); h = sclose(h);
		elseif isstruct(fn{k}) & isfield(fn{k},'EVENT')	
			h = fn{k};
		end	
		eventmatrix = [h.EVENT.POS,h.EVENT.TYP,h.EVENT.CHN,h.EVENT.DUR];
                if k==1,
                        HDR = h; 
                        if length(HDR.FILE)>1,
                                fprintf(2,'Warning Artifact_Selection: does not work correctly for multiple files.\n');
                        end;
                        EVENTMATRIX = eventmatrix; 
		else	
			if isfield(h.EVENT,'SampleRate')
				if isfield(HDR.EVENT,'SampleRate')
					if ~isnan(h.EVENT.SampleRate),
					if HDR.EVENT.SampleRate ~= h.EVENT.SampleRate,
						eventmatrix(:,[1,4]) = eventmatrix(:,[1,4])*HDR.EVENT.SampleRate/h.EVENT.SampleRate;
					end	
					end
				else
					HDR.SampleRate = h.SampleRate;
				end;	 
			end; 
			% merge
			EVENTMATRIX = [EVENTMATRIX; eventmatrix]; 
		end;	
	end;	
end;
%[tmp,ix]=sort(EVENTMATRIX(:,1));
EVENTMATRIX = unique(EVENTMATRIX,'rows');  %$ remove double entries
 
HDR.EVENT.POS = EVENTMATRIX(:,1);
HDR.EVENT.TYP = EVENTMATRIX(:,2);
HDR.EVENT.CHN = EVENTMATRIX(:,3);
HDR.EVENT.DUR = EVENTMATRIX(:,4);

% prepare trigger information 
if ~isfield(HDR,'TRIG')
	HDR.TRIG = HDR.EVENT.POS(HDR.EVENT.TYP==hex2dec('0300'));
end;
HDR.TRIG = sort(HDR.TRIG);
SEL = logical(zeros(length(HDR.TRIG),1));

% define interval
if nargin<2, 
	if HDR.FLAG.TRIGGERED, 
		ti = [1, HDR.SPR];
	else
		ti = [1, max(diff(HDR.TRIG))];
	end;	
else
	if prod(size(t1))==1,
	        ti = [t1,t2];
        elseif prod(size(t1))==2,
	        ti = t1;
        else
                error('invalid time interval');
                return; 
	end;    
	ti = ti*HDR.SampleRate;
end;

if min(diff(HDR.TRIG))<(ti(2)-ti(1))
	fprintf(2,'Warning: trials do overlap.\n');
end;

% prepare artifact information 
ix  = find(bitand(HDR.EVENT.TYP,hex2dec('FFF0'))==hex2dec('0100'));
A.EVENT.POS = [HDR.EVENT.POS(ix); HDR.EVENT.POS(ix) + HDR.EVENT.DUR(ix); inf];   % onset and offset 
A.EVENT.TYP = [ones(length(ix),1); -ones(length(ix),1); 0];			% onset = +1, offset = -1;
[A.EVENT.POS, ix2] = sort(A.EVENT.POS);		%  sort the positions
A.EVENT.TYP = A.EVENT.TYP(ix2);		

% check each trial for an artifact. 
TRIG = [HDR.TRIG(:);inf];
ix1 = 1; ix2 = 1; a = 0; k=0;
P1 = HDR.TRIG(ix1)+ti(1);
P2 = A.EVENT.POS(ix2);
P3 = TRIG(ix1)+ti(2);
while (ix1<=length(HDR.TRIG)) & (ix2<length(A.EVENT.POS))
	k = k+1;
	%[P1,P2,P3]
	if P1<=P2, 
		SEL(ix1) = (P3>P2) | a;
		%fprintf(1,'%6i\t',-1,k,a, ix1,ix2,P1,P2,P3,SEL(ix1));
		ix1 = ix1+1;
		P1  = TRIG(ix1)+ti(1);
		P3  = TRIG(ix1)+ti(2);
	elseif P2<P1, 
		a   = a + A.EVENT.TYP(ix2); 
		%fprintf(1,'%6i\t',-2,k,a,ix1,ix2,P1,P2,P3,SEL(ix1));
		ix2 = ix2+1;
		P2  = A.EVENT.POS(ix2);
	end;
	%fprintf(1,'\n');
end; 

HDR.ArtifactSelection = SEL; 
