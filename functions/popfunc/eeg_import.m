% EEG_IMPORT - Import data files in a variety of supported format
%
% Usage:  
%  >> EEG = eeg_import(fileIn);
%  >> EEG = eeg_import(fileIn, 'key', 'val', ...); 
%
% Inputs:
%   fileIn - [string] input file name. The function used to import the file is 
%            based on the file extension.
%   '.set'  - POP_LOADSET (EEGLAB data format)
%   '.edf'  - POP_BIOSIG (European Data Format)
%   '.bdf'  - POP_BIOSIG (BIOSEMI)
%   '.vhdr' - POP_LOADBV (Brain Vision Exchange Format)
%   '.eeg'  - POP_LOADBV (Brain Vision Exchange Format)
%   '.cnt'  - POP_LOADCNT (Neuroscan)
%   '.mff'  - POP_MFFIMPORT (EGI MFF data format)
%   '.raw'  - POP_READEGI (EGI Binary Simple format)
%   '.nwb'  - POP_NWBIMPORT (Neurodata Without Borders)
%   '.mefd' - POP_MEF3 (MEF iEEG data format)
%   '.ds'   - POP_FILEIO or POP_CTF_READ (see below; CTF MEG data format)
%   '.fif'  - POP_FILEIO (FIF MEG data format)
%   '.fif.gz' - POP_FILEIO (FIF MEG data format)
%
% Optional Inputs:
%  'noevents' - ['on'|'off'] do not save event files. Default is 'off'.
%  'modality' - ['eeg'|'ieeg'|'meg'|'auto'] type of data. 'auto' means the
%               format is determined from the file extension. Default is 'auto'.
%  'ctffunc'  - ['fileio'|'ctfimport'] function to use to import CTF data. 
%               Some data is better imported using the POP_CTF_READ
%               function and some data using POP_FILEIO. Default 'fileio'.
%  'importfunc' - [function handle or string] function to import raw data to
%              EEG format. Default: auto.
%
% Outputs:
%   EEG  - EEGLAB data structure 
% 
% Authors: Arnaud Delorme, SCCN, INC, UCSD, 2024

% Copyright (C) Arnaud Delorme, 2024
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [EEG, modality] = eeg_import(fileIn, varargin)

if nargin < 1
    help eeg_import
    return
end

opt = finputcheck(varargin, {
    'ctffunc'      'string'    { 'fileio' 'ctfimport' }    'fileio'; ...
    'importfunc'   ''  {}    '';
    'importfunc'   ''  {}    '';
    'modality'     'string'  {'ieeg' 'meg' 'eeg' 'auto'}       'auto';
    'noevents'     'string'  {'on' 'off'}    'off' }, 'eeg_import');
if isstr(opt), error(opt); end

[~,~,ext] = fileparts(fileIn);
ext = lower(ext);
if ~isempty(opt.importfunc)
    if strcmpi(opt.modality, 'auto'), opt.modality = 'eeg'; end
    EEG = feval(opt.importfunc, fileIn);
elseif strcmpi(ext, '.bdf') || strcmpi(ext, '.edf')
    if strcmpi(opt.modality, 'auto'), opt.modality = 'eeg'; end
    EEG = pop_biosig(fileIn);
elseif strcmpi(ext, '.vhdr')
    if strcmpi(opt.modality, 'auto'), opt.modality = 'eeg'; end
    [fpathin, fname, ext] = fileparts(fileIn);
    EEG = pop_loadbv(fpathin, [fname ext]);
elseif strcmpi(ext, '.set')
    if strcmpi(opt.modality, 'auto'), opt.modality = 'eeg'; end
    EEG = pop_loadset(fileIn);
elseif strcmpi(ext, '.cnt')
    if strcmpi(opt.modality, 'auto'), opt.modality = 'eeg'; end
    EEG = pop_loadcnt(fileIn, 'dataformat', 'auto');
    datFile = [fileIn(1:end-4) '.dat'];
    if exist(datFile,'file')
        EEG = pop_importevent(EEG, 'indices',1:length(EEG.event), 'append','no', 'event', datFile,...
            'fields',{'DatTrial','DatResp','DatType','DatCorrect','DatLatency'},'skipline',20,'timeunit',NaN,'align',0);
    end
elseif strcmpi(ext, '.mff')
    if strcmpi(opt.modality, 'auto'), opt.modality = 'eeg'; end
    EEG = pop_mffimport(fileIn,{'code'});
elseif strcmpi(ext, '.raw')
    if strcmpi(opt.modality, 'auto'), opt.modality = 'eeg'; end
    EEG = pop_readegi(fileIn);
elseif strcmpi(ext, '.eeg')
    if strcmpi(opt.modality, 'auto'), opt.modality = 'eeg'; end
    [tmpPath,tmpFileName,~] = fileparts(fileIn);
    if exist(fullfile(tmpPath, [tmpFileName '.vhdr']), 'file')
        EEG = pop_loadbv( tmpPath, [tmpFileName '.vhdr'] );
    else
        error('.eeg files not from BrainVision are currently not supported')
    end
elseif strcmpi(ext, '.nwb')
    if strcmpi(opt.modality, 'auto'), opt.modality = 'ieeg'; end
    if ~exist('pop_nwbimport', 'file')
        error('NWB-io plugin not present, please install the plugin first')
    end
    EEG = pop_nwbimport(fileIn, 'importspikes', 'on', 'typefield', 1);
elseif strcmpi(ext, '.mefd')
    if strcmpi(opt.modality, 'auto'), opt.modality = 'ieeg'; end
    modality = 'ieeg';
    if ~exist('pop_MEF3', 'file')
        error('MEF plugin not present, please install the MEF3 plugin first')
    end
    EEG = pop_MEF3(eegFileRaw); % MEF folder
elseif strcmpi(ext, '.fif')
    if strcmpi(opt.modality, 'auto'), opt.modality = 'meg'; end
    EEG = pop_fileio(eegFileRaw); % fif folder
elseif strcmpi(ext, '.gz')
    if strcmpi(opt.modality, 'auto'), opt.modality = 'meg'; end
    gunzip(eegFileRaw);
    EEG = pop_fileio(eegFileRaw(1:end-3)); % fif folder
elseif strcmpi(ext, '.ds')
    if strcmpi(opt.modality, 'auto'), opt.modality = 'meg'; end
    if strcmpi(opt.ctffunc, 'fileio')
        EEG = pop_fileio(eegFileRaw);
    else
        EEG = pop_ctf_read(eegFileRaw);
    end
else
    error('Data format not supported');
end
if strcmpi(opt.noevents, 'on')
    EEG.event = [];
end
modality = opt.modality;