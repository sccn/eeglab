function [signal,H] = tload(FILENAME,TI1,CHAN,EVENTFILE,TI2)
% TLOAD loads and triggers signal data.  
%
% [signal,HDR] = tload(FILENAME, TI, [CHAN,] EVENTFILE, AI)
%
% S = reshape(signal,HDR.size) returns the corresponding 3-dim Matrix 
%
% FILENAME  name of file, or list of filenames, wildcards '*' are supported. 
%	    The files must contain the trigger information. 
% TI	    trigger interval [t1,t2,t3] in seconds, relative to TRIGGER point
%	    The interval [t1,t2] defines the trigger segment, t3 is 
%	    optional and determines the number of NaN's after each trial. 
% CHAN      list of selected channels
%           default=0: loads all channels
% EVENTFILE file of artifact scoring 
% AI	    Artifactinterval [t1,t2] in seconds, relative to TRIGGER point  	
%	    Trials with artifacts within this segment are removed. 
%	    By default AI=TI, AI enables to select the critical period. 
%
%
% see also: SLOAD, SVIEW, SOPEN, ARTIFACT_SELECTION


%	$Id: tload.m,v 1.1 2009-07-07 02:23:48 arno Exp $
%	Copyright (C) 2004-2005 by Alois Schloegl <a.schloegl@ieee.org>
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


if nargin<3; CHAN=0; end;
if isempty(CHAN), CHAN = 0; end;
if (CHAN<1) | ~isfinite(CHAN),
        CHAN=0;
end;

[s,HDR] = sload(FILENAME,CHAN); 
if nargin>3,
        if nargin<5, TI2 = TI1(1:2); end; 
        HDR = artifact_selection({HDR,EVENTFILE},TI2);
end;

H = HDR;
TRIG = HDR.TRIG;
if isfield(HDR,'ArtifactSelection'),
        fprintf(1,'TLOAD: Due to Artifact_Selection, %i (out of %i) trials have been removed from %s \n', sum(HDR.ArtifactSelection),length(TRIG),HDR.FileName);
        TRIG = TRIG(~HDR.ArtifactSelection);
        if isfield(HDR,'Classlabel');
                22,
                H.Classlabel = HDR.Classlabel(~HDR.ArtifactSelection);
        end;
        H.ArtifactSelection = zeros(size(TRIG));
end;

TI1 = TI1*HDR.SampleRate;
if length(TI1)<3, TI1(3)=0; end; 
if HDR.FLAG.TRIGGERED & (any(TI1<1) | any(TI1>HDR.SPR))
	fprintf(2,'Warning TLOAD: data is already triggered - invalid trigger interval\n');
	[signal,sz] = trigg(s,TRIG,1,HDR.SPR,TI1(2)-TI1(1)-HDR.SPR+TI1(3));
	signal = [signal(:,1+end+TI1(1):end),signal(:,1:end+TI1(1))];
else
        [signal,sz] = trigg(s,TRIG,TI1(1),TI1(2),TI1(3));
end;		
signal = signal'; 

H.EVENT = [];
H.size = sz([2,3,1]);
H.FLAG.TRIGGERED = 1; 
H.SPR  = sz(2);
H.DUR  = diff(TI1)/H.SampleRate;
H.NRec = sz(3); 
H.TRIG = (0:H.NRec-1)'*H.SPR; 
        