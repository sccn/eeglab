function [POS,HDR] = stell(HDR)
% STELL returns file position of signal data files
% HDR = stell(HDR)
% returns the location of the HDR-signal file position indicator in the specified file.  
% Position is indicated in Blocks from the beginning of the file.  If -1 is returned, 
% it indicates that the query was unsuccessful; 
% HDR-struct is a struct obtained by SOPEN.
%
% HDR.FILE.POS contains the position of the HDR-Identifier in Blocks
%
% See also: SOPEN, SREAD, SWRITE, SCLOSE, SSEEK, SREWIND, STELL, SEOF


%	$Id: stell.m,v 1.1 2009-07-07 02:23:48 arno Exp $
%	(C) 1997-2005,2007 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/


POS = ftell(HDR.FILE.FID);
if POS<0,
        [HDR.ERROR,HDR.ErrNo] = ferror(HDR.FILE.FID);
        return; 
end;
POS = (POS-HDR.HeadLen)/HDR.AS.bpb;
if (POS<0) || (POS>HDR.NRec*HDR.SPR)
        fprintf(2,'Warning STELL: %s File - invalid file position %i  %i\n', HDR.FileName, POS, HDR.FILE.POS);
elseif ~isfield(HDR.FILE,'POS')
    	fprintf(2,'Error STELL: format %s not supported',HDR.TYPE);
elseif HDR.FILE.POS~=POS,
        fprintf(2,'Warning STELL: %s File position error  %i  %i\n', HDR.FileName, POS, HDR.FILE.POS);
else
	HDR.FILE.POS=POS;	
end;        
 
return; 

% the code below is obsolete 

if strmatch(HDR.TYPE,{'CTF'}),
	POS = (POS-HDR.HeadLen)/HDR.AS.bpb;
	HDR.ERROR = [];
	HDR.ErrNo = 0;

	if (HDR.AS.startrec+HDR.AS.numrec)~=POS,
        	fprintf(2,'Warning STELL: File postion error in EDF/GDF/BDF-toolbox.\n')
                HDR.AS.startrec = POS;
        end;

elseif strmatch(HDR.TYPE,{'AINF','BKR','ISHNE','CNT','EEG','AVG','MIT','RG64','LABVIEW','Nicolet','EGI','SMA','SND','WAV','AIF','CFWB','DEMG','alpha','BCI2000','ET-MEG','Sigma'}),
	POS = (POS-HDR.HeadLen)/HDR.AS.bpb;

elseif strmatch(HDR.TYPE,{'ACQ','AINF','EDF','BDF','EPL','GDF','RDF','SIGIF','BrainVision','EEProbe-CNT','EEProbe-AVR','FIF','native','MFER','TMS32','WG1'}),
	POS = HDR.FILE.POS;
	
else
    	fprintf(2,'Error STELL: format %s not supported',HDR.TYPE);    
end;        

if HDR.FILE.POS~=POS,
        fprintf(2,'Warning STELL: %s File position error  %i  %i\n', HDR.FileName, POS, HDR.FILE.POS);
end;        
HDR.FILE.POS=POS;	
