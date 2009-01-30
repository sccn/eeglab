function h0=get_bbci_regress_eog(fn)
% GET_BBCI_REGRESS_EOG tries to obtain the regression coefficients
%    for EOG correction. 
% 
% Warning: this is a utility function used by biosig. do not use it 
% directly, unless you know what you are doing. At least you are warned. 
%
% See also: SLOAD, IDENTIFY_EOG_CHANNELS, BV2BIOSIG_EVENTS 

%	$Id: get_bbci_regress_eog.m,v 1.1 2009-01-30 06:04:43 arno Exp $
%	Copyright (C) 2006 by Alois Schloegl 
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


	% X-tract information from BBCI data for EOG correction
	if nargin<2,
		CHAN = 0 ; 
	end; 	

	[p,f,e]=fileparts(fn); 

	f0 = dir(fullfile(p,'arte*.vhdr')); 
	if length(f0)<1
		fprintf('error: no artefact file found in %s\n',p); 
		return; 
	elseif length(f0)>1
		fprintf('error: more than one artefact file found in %s\n',p); 
		{fn.name}
		return; 
	end; 

	% load EOG artefacts data
        [s00,h0] = sload(fullfile(p,f0.name));
        
        % define EOG channels
        eogchan = identify_eog_channels(h0); 
        % remove EOG channels for list of corrected channels
	chan = find(~any(eogchan,2)); 
                
        %%%% convert BV (BBCI) events into BioSig Event-Codes
        h0 = bv2biosig_events(h0);         

        % find eye movements in BBCI recordings
        %ix1 = min([strmatch('Augen',h0.EVENT.Desc);strmatch('augen',h0.EVENT.Desc)]);
        ix1 = min(find(bitand(hex2dec('7ff0'),h0.EVENT.TYP)==hex2dec('0430'))); % first eye movement
        %ix2 = min(strmatch('blinzeln',h0.EVENT.Desc)+1); 
        ix2 = max(find(bitand(2^15-1,h0.EVENT.TYP)==hex2dec('0439'))); % end of blinks
	                
        % regression analysis - compute correction coefficients. 
        h0.REGRESS = regress_eog(s00(h0.EVENT.POS(ix1):h0.EVENT.POS(ix2),:),chan,eogchan); 
	
