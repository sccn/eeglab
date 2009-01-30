function [backpanel, frontpanel] = read_ctf_trigger(dataset)

% READ_CTF_TRIGGER reads the STIM channel from a dataset and detects
% the trigger moments and values
%
% [backpanel, frontpanel] = read_ctf_trigger(dataset)
% 
% This returns all samples of the STIM channel, converted to backpanel
% and frontpanel trigger values. Triggers are placed at the rising flank
% of the STIM channel.
%
% Triggers should be at least 9 samples long and are should not overlap 
% each other.
%
% See also READ_CTF_MEG4, READ_CTF_RES4

% Copyright (C) 2003, Robert Oostenveld
%
% $Log: not supported by cvs2svn $
% Revision 1.1  2005/12/06 06:24:23  psdlw
% Alternative functions from the FieldTrip package, which is now released under GPL (so I assume these functions can be committed to the sourceforge cvs)
%
% Revision 1.4  2003/10/14 12:37:17  roberto
% made teh trigshift adaptive to sampling frequency
%
% Revision 1.3  2003/09/12 09:27:11  roberto
% added comment to help about how triggers are supposed to behave
%
% Revision 1.2  2003/05/21 10:59:52  roberto
% fixed bugs in front/backpanel
%
% Revision 1.1  2003/05/19 14:40:32  roberto
% new implementation, replaces inline code in framework/preprocessing
%

[path, file, ext] = fileparts(dataset);
datafile   = fullfile(dataset, [file '.meg4']);
headerfile = fullfile(dataset, [file '.res4']);

% read the header from the raw CTF data
hdr = read_ctf_res4(headerfile);

% number of samples to shift the assesment of the trigger value
% this is needed because it takes some time for the rising flank to get to the correct value
trigshift = fix(hdr.Fs * 9/1200);

% read the stimulus channel from raw CTF data
stimindx = find(strcmp(hdr.label, 'STIM'));
stim =  read_ctf_meg4(datafile, hdr, 1, hdr.nTrials*hdr.nSamples, stimindx);

% determine the precise timing of the triggers
upflank = [0 (diff(stim)>0 & stim(1:(end-1))==0)];
trigger = upflank(1:(end-trigshift)).*stim((1+trigshift):end);

% determine the triggers on the backpanel, only take highest 16 bits
backpanel = fix(trigger / 2^16);

% determine the triggers on the frontpanel, only take lowest 16 bits
frontpanel = double(bitand(uint32(trigger), 2^16-1));

