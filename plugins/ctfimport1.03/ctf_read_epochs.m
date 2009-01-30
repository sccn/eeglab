function [data,ctf] = ctf_read_epochs(folder,varargin);

% ctf_read_epochs - extract epochs from a ctf dataset (.ds)
%
% [data,ctf] = ctf_read_epochs(folder,varargin)
%
% This function will extract epochs for specific markers.  The inputs are
% given as 'parameter' => 'value' pairs, where the following parameters and
% values are recognised:
% 'channels' => 'eeg','meg','ref','vc','other','all',[a:b]
% 'markers' => {'marker names'}
% 'window' => [start_sec,end_sec]
% 'trials' => 'all',[a:b]
%
%      <>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%      <                                                       >
%      <                      DISCLAIMER:                      >
%      <                                                       >
%      <  THIS PROGRAM IS INTENDED FOR RESEARCH PURPOSES ONLY. >
%      <  THIS PROGRAM IS IN NO WAY INTENDED FOR CLINICAL OR   >
%      <                     OFFICIAL USE.                     >
%      <                                                       >
%      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>
%

% $Revision: 1.1 $ $Date: 2009-01-30 03:49:27 $

% Copyright (C) 2004  Darren L. Weber
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% Modified: 01/2004, Darren.Weber_at_radiology.ucsf.edu
%                    - modified from NIH code readepochs.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.1 $';
fprintf('\nCTF_READ_EPOCHS [v %s]\n',ver(11:15)); tic;

if ~exist('folder','var'),
  ctf = ctf_folder;
else
  ctf = ctf_folder(folder);
end

ctf = ctf_read_res4(ctf.folder);

% defaults
CHAN = ctf.sensor.index.meg;
TIME = 'alltimes';
TRIALS = [1:ctf.setup.number_trials];
markers = [];

for i = 1:size(varargin,2)
    if iscell(varargin{i})
        if strmatch(varargin{i},'markers'),
            markers = varargin{i+1};
        end    
    elseif ischar(varargin{i}),
        switch lower(varargin{i}),
            case 'channels',
                if isnumeric(varargin{i+1}),
                    CHAN = varargin{i+1};
                else
                    switch lower(varargin{i+1}),
                        case 'eeg',
                            CHAN = ctf.sensor.index.eeg;
                        case 'meg',
                            CHAN = ctf.sensor.index.meg;
                        case 'ref',
                            CHAN = ctf.sensor.index.ref;
                        case 'other',
                            CHAN = ctf.sensor.index.other;
                        case 'vc',
                            CHAN = ctf.sensor.index.vc;
                        case 'all',
                            CHAN = [1:ctf.setup.number_channels];
                    end
                end
            case 'window',
                if isnumeric(varargin{i+1}),
                    TIME = varargin{i+1};
                else
                    error('window argument must be numeric array');
                end
            case 'trials',
                if isnumeric(varargin{i+1}),
                    TRIALS = varargin{i+1};
                else
                    error('trials argument must be numeric array');
                end
        end
    end
end

if isempty(markers),
    
    % If not markers are specified, we just rearrange the ctf.data cell
    % array into a 3D matrix
    
    ctf = ctf_read_meg4(folder,ctf,CHAN,TIME,TRIALS);
    
    for i = 1:size(ctf.data,1),
        
        epochs{1}(:,:,i) = ctf.data{i};
        
    end
    
else
    
    % read the marker file
    ctf = ctf_read_markerfile(ctf.folder,ctf);
    
    if isempty(ctf.markers),
        error('No markers');
    end
    
    number_markers = length(markers);
    
    % initialise epoch cell array to hold an ordered set of epochs for each
    % marker selected
    epochs = cell(number_markers,1);
    
    for marker = 1:number_markers,
        
        mk = ismember(ctf.markers.marker_names,markers(marker));
        
        keyboard
        
        ctf.markers.marker_names(mk);
        
        nsamp = ctf.markers.number_samples(mk)
        
        nss=0;
        
        for ns=1:nsamp
            
            tr = ctf.markers.trial_times{mk}(ns,1);
            
            if ismember(tr,TRIALS)
                
                nss=nss+1;
                
                tim=ctf.markers.trial_times{mk}(ns,2);
                
                times=TIME+tim;
                
                ctf = ctf_read_meg4(folder,ctf,CHAN,TIME,tr);
                
                temp=ctf.data{1};
                
                epochs{marker}(:,:,nss)=temp;
                
            end
        end
    end
end

data = struct(...
    'folder',ctf.folder,...
    'setup', ctf.setup,...
    'sensor',ctf.sensor,...
    'epochs',epochs );

return
