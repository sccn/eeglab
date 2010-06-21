% pop_loadbva() - import a Matlab file from brain vision analyser
%                  software.
%
% Usage:
%   >> OUTEEG = pop_loadbva( filename );
%
% Inputs:
%   filename       - file name
%
% Outputs:
%   OUTEEG         - EEGLAB data structure
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, Dec 2003
%
% See also: eeglab()

% Copyright (C) 2003 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [EEG, com] =  pop_loadbva(filename)
    
    EEG = [];
    com = '';
    
    if nargin < 1
        [tmpfilename, filepath] = uigetfile('*.mat;*.MAT', 'Choose a Matlab file from Brain Vision Analyser -- pop_loadbva');
        if tmpfilename == 0 return; end;
        filename = [ filepath tmpfilename ];
    end;
    
    disp('Importing data');
    EEG = eeg_emptyset;
    bva = load(filename, '-mat');
    allfields = fieldnames(bva);
    chanstruct = bva.Channels;
    channames  = lower({ chanstruct.Name });
    for index = 1:length(allfields)
        switch lower(allfields{index})
         case { 't' 'markercount' 'markers' 'samplerate' 'segmentcount' 'channelcount' 'channels' },
         otherwise
          count1 = strmatch(lower(allfields{index}), channames, 'exact');
          count2 = strmatch(lower(allfields{index}(2:end)), channames, 'exact');
          if ~isempty(count1) | ~isempty(count2)
              count = [ count1 count2 ];
              count = count(1);
          else
              disp(['Warning: channel ''' lower(allfields{index}) ''' not in channel location structure']);
              count = length(chanstruct)+1;
              chanstruct(end+1).Name   = allfields{index};
              chanstruct(end+1).Phi    = [];
              chanstruct(end+1).Theta  = [];
              chanstruct(end+1).Radius = [];
          end;
          if bva.SegmentCount > 1
              EEG.data(count,:,:) = getfield(bva, allfields{index})';
              bva = rmfield(bva, allfields{index});
          else
              EEG.data(count,:) = getfield(bva, allfields{index});
              bva = rmfield(bva, allfields{index});
          end;
        end;
    end;
    
    EEG.nbchan   = size(EEG.data,1);
    EEG.srate    = bva.SampleRate;
    EEG.xmin     = bva.t(1)/1000;
    EEG.xmax     = bva.t(end)/1000;
    EEG.pnts     = size(EEG.data,2);
    EEG.trials   = size(EEG.data,3);
    EEG.setname  = 'Brain Vision Analyzer file';
    EEG.comments = [ 'Original file: ' filename ];

    % convert channel location structure
    % ----------------------------------
    disp('Importing channel location information');
    for index = 1:length(chanstruct)
        EEG.chanlocs(index).labels = chanstruct(index).Name;
        if chanstruct(index).Radius ~= 0
            EEG.chanlocs(index).sph_theta_besa = chanstruct(index).Theta;
            EEG.chanlocs(index).sph_phi_besa   = chanstruct(index).Phi;
            EEG.chanlocs(index).sph_radius     = chanstruct(index).Radius;
        else
            EEG.chanlocs(index).sph_theta_besa = [];
            EEG.chanlocs(index).sph_phi_besa   = [];
            EEG.chanlocs(index).sph_radius     = [];
        end;
    end;
    EEG.chanlocs = convertlocs(EEG.chanlocs, 'sphbesa2all');
    EEG.chanlocs = rmfield(EEG.chanlocs, 'sph_theta_besa');
    EEG.chanlocs = rmfield(EEG.chanlocs, 'sph_phi_besa');
    
    % convert event information
    % -------------------------
    disp('Importing events');
    index = 0;
    for index1 = 1:size(bva.Markers,1)
        for index2 = 0:size(bva.Markers,2)-1
            if ~isempty(bva.Markers(index2*size(bva.Markers,1)+index1).Description)
                index = index + 1;
                EEG.event(index).type        = bva.Markers(index2*size(bva.Markers,1)+index1).Description;
                EEG.event(index).latency     = bva.Markers(index2*size(bva.Markers,1)+index1).Position;
                EEG.event(index).Points      = bva.Markers(index2*size(bva.Markers,1)+index1).Points;
                try
                    EEG.event(index).bvatype     = bva.Markers(index2*size(bva.Markers,1)+index1).Type;     
                    EEG.event(index).description = bva.Markers(index2*size(bva.Markers,1)+index1).Description;                
                catch, end;
                try
                    EEG.event(index).chan        = bva.Markers(index2*size(bva.Markers,1)+index1).Chan;
                catch, end;
                try
                    EEG.event(index).channelnumber = bva.Markers(index2*size(bva.Markers,1)+index1).ChannelNumber;
                catch, end;
                if bva.SegmentCount > 1
                    EEG.event(index).epoch       = index1;
                    EEG.event(index).latency     = bva.Markers(index2*size(bva.Markers,1)+index1).Position+(index1-1)*EEG.pnts;
                else
                    EEG.event(index).latency     = bva.Markers(index2*size(bva.Markers,1)+index1).Position;
                end;
            end;
        end;
    end;
    EEG = eeg_checkset(EEG, 'makeur');
    EEG = eeg_checkset(EEG, 'eventconsistency');
    
    com = sprintf('EEG = pop_loadbva(''%s'');', filename);
