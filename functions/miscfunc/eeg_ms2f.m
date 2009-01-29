% eeg_ms2f() - convert epoch latency in ms to nearest epoch frame number
%
% Usage:
%         >> outf = eeg_ms2f(EEG,ms);
% Inputs:
%         EEG - EEGLAB data set structure
%         ms  - epoch latency in milliseconds
% Output:
%         outf - nearest epoch frame to the specified epoch latency
%
% Author: Scott Makeig, SCCN/INC, 12/05

function outf = eeg_ms2f(EEG,ms)
        ms = ms/1000;
        if ms < EEG.xmin | ms > EEG.xmax
           error('time out of range');
        end
        outf = 1+round((EEG.pnts-1)*(ms-EEG.xmin)/(EEG.xmax-EEG.xmin));
        % else
        % [tmp outf] = min(abs(EEG.times-ms)); 
