% pop_rejspec() - rejection of artifact in a dataset using 
%               thresholding of frequencies in the data.
%
% Usage:
%   >>  pop_rejspec(INEEG, typerej); % pops-up
%   >> [OUTEEG, Indexes] = pop_rejspec( INEEG, typerej, electrodes, ...
%         negthresh, posthresh, startfreq, endfreq, superpose, reject);
%
% Inputs:
%   INEEG      - input dataset
%   typerej    - type of rejection (0 = independent components; 1 = eeg
%              data). Default is 1. For independent components, before
%              thresholding, the activity is renormalized for each 
%              component.
%   electrodes - [e1 e2 ...] electrodes (number) to take into 
%              consideration for rejection
%   negthresh  - negative threshold limit in mV (can be an array if 
%              several electrodes; if less numbe  of values than number 
%              of electrodes the last value is used for the remaining 
%              electrodes)
%   posthresh  - positive threshold limit in mV (same syntax as negthresh
%   startfreq  - starting frequency in Hz (same syntax  as negthresh)
%   endfreq    - ending frequency in Hz (same syntax  as negthresh).
%              Starfreq and endfreq define the frequncy range for
%              rejection.
%   superpose  - 0=do not superpose pre-labelling with previous
%              pre-labelling (stored in the dataset). 1=consider both
%              pre-labelling (using different colors). Default is 0.
%   reject     - 0=do not reject labelled trials (but still store the 
%              labels. 1=reject labelled trials. Default is 0.
%
% Outputs:
%   OUTEEG     - output dataset with updated spectrograms
%   Indexes    - index of rejected sweeps
%     when eegplot is called, modifications are applied to the current 
%     dataset at the end of the call of eegplot (when the user press the 
%     button 'reject').
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eegtresh(), eeglab(), eegplot(), pop_rejepoch()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: not supported by cvs2svn $
% Revision 1.11  2002/08/12 02:32:24  arno
% inputdlg2
%
% Revision 1.10  2002/08/07 22:38:42  arno
% editing header
%
% Revision 1.9  2002/07/31 01:02:13  arno
% debugging
%
% Revision 1.8  2002/07/30 23:37:51  arno
% debugging
%
% Revision 1.7  2002/07/30 23:33:47  arno
% new rejection type
%
% Revision 1.6  2002/07/30 15:17:52  arno
% debug display
%
% Revision 1.5  2002/07/30 15:17:16  arno
% debugging
%
% Revision 1.4  2002/07/26 18:09:31  arno
% debugging
%
% Revision 1.3  2002/07/26 18:04:48  arno
% debugging
%
% Revision 1.2  2002/07/26 16:44:35  arno
% switching icacomp
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 03-07-02 added srate argument to eegplot call -ad
% 03-08-02 reworked spectrum to save space & add eeglab options -ad

function [EEG, I1, com] = pop_rejspec( EEG, icacomp, elecrange, negthresh, posthresh, ...
   						startfreq, endfreq, superpose, reject, topcommand);

I1 = [];
com = '';
if nargin < 1
   help pop_rejspec;
   return;
end;  
if nargin < 2
   icacomp = 1;
end;  

if icacomp == 0
	if isempty( EEG.icasphere )
	    ButtonName=questdlg( 'Do you want to run ICA now ?', ...
                         'Confirmation', 'NO', 'YES', 'YES');
    	switch ButtonName,
        	case 'NO', disp('Operation cancelled'); return;   
        	case 'YES', [ EEG com ] = pop_runica(EEG);
    	end % switch
	end;
end;	
if nargin < 3

	% which set to save
	% -----------------
	promptstr   = { fastif(icacomp==0, 'Component (number; ex: 2 4 5):', 'Electrode (number; ex: 2 4 5):'), ...
					'Lower limit(s) (dB):', ...
					'Upper limit(s) (dB):', ...
					'Low frequency(s) (Hz):', ...
					'High frequency(s) (Hz):', ...
               		'Display with previous rejection', ...
         			'Actually reject marked trial (YES or NO)' };
	inistr      = { ['1:' int2str(EEG.nbchan)], ...
					'-10', ...
					'10', ...
					'15', ...
					'30', ...
               		'NO', ...
            		'NO' };

	result       = inputdlg2( promptstr, fastif(~icacomp, 'Reject by component spectra -- pop_rejspec()', ...
											   'Reject by data spectra -- pop_rejspec()'), 1,  inistr, 'pop_rejspec');
	size_result  = size( result );
	if size_result(1) == 0 return; end;
	elecrange    = result{1};
	negthresh    = result{2};
	posthresh    = result{3};
	startfreq    = result{4};
	endfreq      = result{5};
	switch lower(result{6}), case 'yes', superpose=1; otherwise, superpose=0; end;
	switch lower(result{7}), case 'yes', reject=1; otherwise, reject=0; end;
end;

if isstr(elecrange) % convert arguments if they are in text format 
	calldisp = 1;
	elecrange = eval( [ '[' elecrange ']' ]  );
	negthresh = eval( [ '[' negthresh ']' ]  );
	posthresh = eval( [ '[' posthresh ']' ]  );
	startfreq = eval( [ '[' startfreq ']' ]  );
	endfreq   = eval( [ '[' endfreq ']' ]  );
else
	calldisp = 0;
end;

sizewin = 2^nextpow2(EEG.pnts);
if icacomp == 1
    [allspec, Irej, rejE, freqs ] = spectrumthresh( EEG.data, EEG.specdata, ...
							elecrange, EEG.srate, negthresh, posthresh, startfreq, endfreq);
else
    % test if ICA was computed
    % ------------------------
    eeg_options; % changed from eeglaboptions 3/30/02 -sm
 	if option_computeica  
    	icaacttmp = EEG.icaact(elecrange, :, :);
	else
        icaacttmp = (EEG.icaweights(elecrange,:)*EEG.icasphere)*reshape(EEG.data, ...
												  EEG.nbchan, EEG.trials*EEG.pnts);
        icaacttmp = reshape( icaacttmp, length(elecrange), EEG.pnts, EEG.trials);
    end;
   	try, oldspec   = EEG.specicaact(elecrange, :, :); catch, oldspec = []; end;
    [allspec, Irej, rejE, freqs ] = spectrumthresh( icaacttmp, oldspec, 1:length(elecrange), ...
								EEG.srate, negthresh, posthresh, startfreq, endfreq);
end;

fprintf('%d channel selected\n', size(elecrange(:), 1));
fprintf('%d/%d trials marked for rejection\n', length(Irej), EEG.trials);
rej = zeros( 1, EEG.trials);
rej(Irej) = 1;

if calldisp
	nbpnts = fastif( icacomp == 1, size(EEG.specdata,2), size(EEG.specicaact,2));
    if icacomp == 1 macrorej  = 'EEG.reject.rejfreq';
        			macrorejE = 'EEG.reject.rejfreqE';
    else			macrorej  = 'EEG.reject.icarejfreq';
        			macrorejE = 'EEG.reject.icarejfreqE';
    end;
	colrej = EEG.reject.rejfreqcol;
	eeg_rejmacro; % script macro for generating command and old rejection arrays
	if icacomp == 1
		eegplot(EEG.data(elecrange,:,:), 'winlength', 5, 'position', [100 800 800 500], ...
				'limits', [EEG.xmin EEG.xmax]*1000, 'xgrid', 'off', 'tag', 'childEEG' );
	else
		eegplot(icaacttmp, 'winlength', 5, 'position', [100 800 800 500], 'limits', ...
				[EEG.xmin EEG.xmax]*1000 , 'xgrid', 'off', 'tag', 'childEEG' );
	end;	
	eegplot( allspec, 'srate', EEG.srate, 'freqlimits', [1 EEG.srate/2], 'command', ...
			 command, 'children', gcf, eegplotoptions{:}); 
end;

% store variables if necessary
% -----------------------------
eeg_options; % changed from eeglaboptions 3/30/02 -sm
if option_keepdataset
    if icacomp == 1, EEG.specdata = zeros(EEG.nbchan, size(allspec,2), size(allspec,3));
                     EEG.specdata(elecrange, :, :) = allspec;
    else,            EEG.specicaact = zeros(EEG.nbchan, size(allspec,2), size(allspec,3));
                     EEG.specicaact(elecrange, :, :) = allspec;
    end;
end;
    
com = [com sprintf('%s = pop_rejspec( %s, %s);', inputname(1), ...
   inputname(1), vararg2str({icacomp, elecrange,  negthresh, posthresh, startfreq, endfreq, superpose, reject })) ]; 
if nargin < 3 & nargout == 2
	I1 = com;
end;

return;

% test if one must recompute the spectrum for some electrodes
% -----------------------------------------------------------
function res = testgoinloop( dat, elec)
    res = 0;
    if isempty(dat), res =1; return; end;
	for index=elec
	   if index > size(dat,1), res = 1; return; end;
	   if all(dat(index,1:10,1) == 0)  res = 1; return; end;
    end;

% compute spectrum and reject artifacts
% -------------------------------------
function [specdata, Irej, Erej, freqs ] = spectrumthresh( data, specdata, elecrange, srate, negthresh, posthresh, startfreq, endfreq);
	% compute the fft if necessary - old version
	%freqs = EEG.srate*[1, sizewin]/sizewin/2;
	%EEG.specdata = fft( EEG.data-repmat(mean(EEG.data,2), [1 EEG.pnts 1]), sizewin, 2);
	%EEG.specdata = EEG.specdata( :, 2:sizewin/2+1, :);
	%EEG.specdata = 20*log2(abs( EEG.specdata ).^2);

	[tmp freqs] = pmtm( data(1,:,1) ); % just to get the frequencies 	
    freqs = [freqs(1) freqs(end)]*srate;	
    
    Irej    = [];
    Erej    = zeros(size(data,1), size(data,3));
	fprintf('Computing spectrum (using slepian tapers; done only once):\n');    
    for index = 1:length(elecrange)
	   if testgoinloop( specdata, index )
	      fprintf('%d ', elecrange(index));    
		  for indextrials = 1:size(data,3)
			[ tmpspec(1,:,indextrials) freqs] = pmtm( data(elecrange(index),:,indextrials) );
		  end;
          freqs = [freqs(1) freqs(end)]*srate;	
 		  tmpspec  = 10*log(tmpspec);
 		  tmpspec  = tmpspec - repmat( mean(tmpspec,3), [1 1 size(data,3)]);
          specdata(index,:,:) = tmpspec;
	   else
	       tmpspec = specdata(index,:,:);
	   end;

       % perform the rejection
       % ---------------------	
	   [I1 Itmprej NS Etmprej] = eegthresh( tmpspec, size(tmpspec,2), 1, negthresh, posthresh, ...
						freqs, startfreq, max(max(freqs), endfreq));
 	   Irej = union(Irej, Itmprej);
 	   Erej(elecrange(index),Itmprej) = Etmprej;
	end;
	fprintf('\n');    
