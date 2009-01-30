function [status]=seof(HDR)
% SEOF checks for end of signal-file
%    status = seof(HDR)
%
% returns 1 if End-of-EDF-File is reached
% returns 0 otherwise
%
% See also: SOPEN, SREAD, SWRITE, SCLOSE, SSEEK, SREWIND, STELL, SEOF

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the  License, or (at your option) any later version.

%	$Id: seof.m,v 1.1 2009-01-30 06:04:42 arno Exp $
%	(C) 1997-2005,2007,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/


%status=feof(HDR.FILE.FID);  % does not work properly
%if HDR.FILE.POS~=HDR.AS.startrec+HDR.AS.numrec;
        
if strmatch(HDR.TYPE,{'CTF','RDF','EEG','AVG','SIGIF'}),
	%status=feof(EDF.FILE.FID);  % does not work properly
	%if EDF.FILE.POS~=EDF.AS.startrec+EDF.AS.numrec;
        status = (HDR.FILE.POS >= HDR.NRec);
	
elseif strmatch(HDR.TYPE,{'RG64','LABVIEW','Nicolet'}),
	status = (HDR.FILE.POS >= (HDR.AS.endpos-HDR.HeadLen));

elseif strmatch(HDR.TYPE,{'ACQ','AINF','BDF','BKR','BrainVision','CNT','CTF','EDF','ET-MEG','GDF','MIT','SMA','CFWB','DEMG','EEProbe-CNT','EEProbe-AVR','MFER','alpha','native','SCP','BCI2000','TMS32','WG1','Sigma'}),
	status = (HDR.FILE.POS >= HDR.SPR*HDR.NRec);

elseif strmatch(HDR.TYPE,{'EGI'}),
        if HDR.FLAG.TRIGGERED,
	        status = (HDR.FILE.POS >= HDR.NRec);
        else        
                status = (HDR.FILE.POS >= HDR.SPR);
        end;

elseif strmatch(HDR.TYPE,{'FIF'}),
        [buf, status] = rawdata('next');
        status = strcmp(status,'eof');
        
else
	status=feof(HDR.FILE.FID);
end;
