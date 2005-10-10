% Usage:    
%   >> [logersp, logfreqs] = cls_readersp(ALLEEG, abset, component);  
%   This functions returns the log scaled ERSP of an ICA component. 
%   The information is loaded from a float file, which a pointer 
%   to is saved in the EEG dataset. The float file was created
%   by the pre - clustering function cls_ersp. 
%
% cls_readersp() - Given the ALLEEG structure, a specific EEG dataset index, 
% and a specific component, the function returns the (frequency) log scaled 
% ERSP of that ICA component. 
% The ERSP of the dataset ICA components is assumed to be saved in a float 
% file, the EEG dataset include a pointer to this file. If such a float file doesn't exist,
% you can use the cls_ersp() function to create it, or use the pre - clustering functions
% that call it: pop_preclust, eeg_preclust & eeg_createdata. The ERSP information  
% is specific to the input variables used to compute it, the frequency range, 
% timewindow, resolution, confidence level and the wavelet type (FFT / wavelet
%  cycles), see timef for details. 
% Along with the ERSP of the selected ICA component the function returns  
% the log frequency vector of the ERSP samples. 
%
%
% Inputs:
%   ALLEEG     - an EEGLAB data structure, which holds EEG sets (can also be one EEG set). 
%                      ALLEEG must contain the dataset of interest (the setind).
%   setind         -  [integer] an index of an EEG dataset in the ALLEEG
%                      structure, for which to get the component log scaled ERSP.
%   component - [integer] a component index in the selected EEG dataset for which 
%                      an ERSP will be returned. 
%
% Outputs:
%   logesrp    - the log scaled ERSP of the requested ICA component in the
%                      selected dataset. This is frequencies (log scaled)
%                      by time, time - frequency decomposition of the ICA
%                      activations.
%   logfreqs   - a vector of the frequencies points in log scale. 
%
%  See also  cls_ersp, pop_preclust, eeg_preclust, eeg_createdata           
%
% Authors:  Hilit Serby, SCCN, INC, UCSD, February, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, October 11, 2004, hilit@sccn.ucsd.edu
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

function [logersp, logfreqs] = cls_readersp(ALLEEG, abset, comp);

for k = 1: length(abset)    
	params = ALLEEG(abset(k)).etc.icaerspparams;
	tlen = length(params.times);
	flen = length(params.freqs);
	try
        erspall{k} =  floatread([ ALLEEG(abset(k)).filepath ALLEEG(abset(k)).etc.icaersp], [flen (tlen+3)],[], flen*(tlen+3)*(comp-1));
        %logersp =  floatread([ ALLEEG(abset).filepath ALLEEG(abset).etc.icalogersp], [flen tlen],[], flen*tlen*(comp-1));
        %erspall =  floatread([ ALLEEG(abset).filepath ALLEEG(abset).etc.icaersp], [flen tlen+2],[], flen*(tlen+3)*(comp-1));
        %ersp = erspall(:, 1:tlen);
        %erspboot = erspall(:,tlen+1:tlen+2).';
	catch
        try
            erspall{k} =  floatread([ ALLEEG(abset(k)).filepath '/' ALLEEG(abset(k)).etc.icaersp], [flen (tlen+3)],[], flen*(tlen+3)*(comp-1));
        catch
            try
                erspall{k} =  floatread([ ALLEEG(abset(k)).filepath '\' ALLEEG(abset(k)).etc.icaersp], [flen (tlen+3)],[], flen*(tlen+3)*(comp-1));
            catch 
                warndlg2(['cls_readersp: file '  ALLEEG(abset(k)).etc.icaersp ' was not found in path ' ALLEEG(abset(k)).filepath], 'Abort - cls_readersp' ); 
                return;
            end
        end
	end
end	
% compute average baseline across conditions
if length(abset) > 1 
	ave_baseline = 0; % mean baseline for requested component across conditions
	for cond = 1:length(abset)
        ave_baseline = ave_baseline + erspall{cond}(:,(tlen+3));
	end
	ave_baseline = ave_baseline./length(abset);
    for cond = 1:length(abset)
        old_base =  erspall{cond}(:,(tlen+3));
        % add back former baseline to the ERSP 
        erspall{cond}(:,1:tlen) = erspall{cond}(:,1:tlen) + 10*log10(repmat(old_base,[1 tlen]));
        % subtruct average baseline from the ERSP
        erspall{cond}(:,1:tlen) = erspall{cond}(:,1:tlen) - 10*log10(repmat(ave_baseline,[1 tlen]));
        % add back former baseline to the ERSPboot 
        erspall{cond}(:,(tlen+1):(tlen+2)) = erspall{cond}(:,(tlen+1):(tlen+2)) + 10*log10(repmat(old_base,[1 2]));
        % subtruct average baseline from the ERSPboot
        erspall{cond}(:,(tlen+1):(tlen+2)) = erspall{cond}(:,(tlen+1):(tlen+2)) - 10*log10(repmat(ave_baseline,[1 2]));  
    end
end
if ~isempty(params.timewindow)
    if params.timewindow(1) > params.times(1) | params.timewindow(end) < params.times(end)
        maxind = max(find(params.times <= params.timewindow(end)));
        minind = min(find(params.times >= params.timewindow(1)));
    else
        minind = 1;
        maxind = tlen;
    end
else
    minind = 1;
    maxind = tlen;
end
if ~isempty(params.freqrange)
    if params.freqrange(1) > exp(1)^params.logfreqs(1) | params.freqrange(end) < exp(1)^params.logfreqs(end)
        fmaxind = max(find(exp(1).^(params.logfreqs) <= params.freqrange(end)));
        fminind = min(find(exp(1).^(params.logfreqs) >= params.freqrange(1)));
    else
        fminind = 1;
        fmaxind = flen;
    end
else
    fminind = 1;
    fmaxind = flen;
end


% Mask ERSP
for cond  = 1:length(abset)
    erspboot = erspall{cond}(:,(tlen+1):(tlen+2)); 
    ersp = erspall{cond}(:,1:tlen); 
	minersp= repmat(erspboot(:,1),1,size(ersp,2));
	maxersp= repmat(erspboot(:,2),1,size(ersp,2));
	ersp(find(ersp<maxersp & ersp>minersp))  = 0;
    ersp = ersp(fminind:fmaxind,minind:maxind);
    logersp(:,:,cond) = ersp;
end

logfreqs = params.logfreqs(fminind:fmaxind);


