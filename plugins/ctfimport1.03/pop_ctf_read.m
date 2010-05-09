% pop_ctf_read() - read CTF file as EEGLAB dataset
%
% Usage:
%   >> EEGOUT = pop_ctf_read; % pop up graphic interface
%   >> EEGOUT = pop_ctf_read(folder);
%   >> EEGOUT = pop_ctf_read(folder, chans, time, trials);
%
% Inputs:
%   folder   - [string]  EEGLAB figure
%   chans    - [integer array or string] see ctf_read()
%   time     - [float array or string] see ctf_read()
%   trials   - [integer array or string] see ctf_read()
%
%
% Author: Arnaud Delorme (SCCN, UCSD) and Daren Weber (DNL, UCSF)
%
% See also: ctf_read(), ctf_readmarkerfile()

% Copyright (C) 2003 Arnaud Delorme, SCCN, UCSD, arno@salk.edu
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

% from ctf2eeglab - script to convert and save ctf .ds into eeglab .set data
%
% The script uses a ctf struct in the matlab workspace or a GUI prompt to
% load a CTF .ds folder, then converts the ctf data into an EEGLAB EEG
% struct, saving the resulting dataset into an EEEGLAB .set file, located
% in the same path as the ctf .ds folder.  The GUI prompt for the CTF .ds
% folder also provides access to definition of the channels, time and
% trials to load.
%
% Licence:  GNU GPL, no express or implied warranties
% Modified: 01/2004, Darren.Weber_at_radiology.ucsf.edu
%                    - developed in collaboration with Fredrick Carver of
%                    the NIH, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [EEG, com] = pop_ctf_read(orifolder, varargin)

  EEG = [];
  com = '';
  if nargin < 1
      ctf = ctf_folder;
      %tmp = ctf_folder(tmp, []);
      %tmp = tmp.folder;
      %tmp = uigetfile('*.*','Select .res4 file');
      %if tmp == 0 return; end;
      
      % read info and prompt
      % --------------------
      listchan = { 'all' 'eeg' 'meg' 'ref' 'other' 'megeeg'};
      disp('Reading file info...');
      ctf = ctf_read_res4(ctf.folder, 0);
      uigeom       = { [2 1] [2 1] [2 1] [2 1] };
      uilist       = { { 'style' 'text' 'string' 'Channels group subset' } ...
                     { 'style' 'list' 'string' strvcat(listchan{:}) } ...
                     { 'style' 'text' 'string' 'or channel indices (overwrite above)' } ...
                     { 'style' 'edit' 'string' '' } ...
                     { 'style' 'text' 'string' 'Time range (default all)' } ...
                     { 'style' 'edit' 'string' '' } ...
                     { 'style' 'text' 'string' sprintf('Trial range (default all [1:%d])', ctf.setup.number_trials) } ...
                     { 'style' 'edit' 'string' '' } };
      result = inputgui( uigeom, uilist, 'pophelp(''pop_ctf_read'')', 'Load a CTF dataset');
      if length( result ) == 0 return; end;
      
      % decode inputs
      % -------------
      options{1} = [1:ctf.setup.number_channels];
      if result{1} > 1
          options{1} = listchan{result{1}};
      end;
      if ~isempty(result{2})
          options{1} = eval( [ '[' result{2} ']' ] );
      end;
      if ~isempty(result{3})
          options{2} = eval( [ '[' result{3} ']' ] );
      end;
      if ~isempty(result{4})
          options{3} = eval( [ '[' result{4} ']' ] );
      end;
  else 
      options = varargin;
  end;
  
  alltrials = [];
  if length(options) == 3
      alltrials = options{3};
  end;
  if length(options) > 1 & isempty(options{2})
      options{2} = 'all';
  end;
  
  % read the data
  % -------------
  ctf = ctf_read(ctf.folder,options{:});
  
  % check if the data is averaged
  % -----------------------------
  if ctf.setup.number_trials_averaged > 0,
      warning('this .ds folder is averaged, removing stdev');
  end

  % check if the data is greater than 500 Mb
  data_size = ctf.setup.number_samples * ctf.setup.number_channels * ctf.setup.number_trials;
  data_bytes = data_size * 8;
  if data_bytes > 5e9, warning('data is greater than 500 Mb'); end
  clear data_size data_bytes;

  % ctf.data is a 3D matrix of samples(time) x channels x trials
  % whereas EEGLAB has a 3D data matrix with channels X samples X trials
  if ndims(ctf.data) < 3
      data = ctf.data';
  else
      try, 
          data = permute(ctf.data, [2 1 3]);
      catch,
          disp('Not enough memory, trying to transpose data by saving it');
          fid = fopen('tmpmeg.raw', 'w')
          if fid == -1,
              error('Cannot create file in current folder');
          end;
          for index = i:size(ctf.data,2)
              fwrite(fid,squeeze(ctf.data(:,index,:)),'float');
          end; 
          fclose(fid);
          ctf.data = [];
          EEG.data = floatread('tmpmeg.raw');
          delete('tmpmeg.raw');
      end;
  end;
  
  % import the data into the EEGLAB EEG struct
  
  [DSpath,DSfile,DSext] = fileparts(ctf.folder);
  
  EEG = eeg_emptyset;
  EEG.setname = DSfile;
  
  % ---
  % These fields now contain the name of the dataset *once*
  % it has been saved (so reamin empty before the dataset 
  % has been saved).
  %EEG.filename = [DSfile,'.set'];
  %EEG.filepath = DSpath;
  % ---

  EEG.comments = [ 'Original folder: ' ctf.folder ];
  %EEG.comments = ctf.setup.run_description';
  EEG.pnts = ctf.setup.number_samples;
  EEG.nbchan = ctf.setup.number_channels;
  EEG.trials = ctf.setup.number_trials;
  EEG.srate = ctf.setup.sample_rate;
  EEG.xmin = ctf.setup.start_sec;
  EEG.xmax = ctf.setup.end_sec;
  EEG.data = data;
  EEG.ref = 'common';

  try,
      for index=1:ctf.setup.number_channels,
          EEG.chanlocs(index).labels = ctf.sensor.label{index};
          EEG.chanlocs(index).X      = ctf.sensor.location(1,index);
          EEG.chanlocs(index).Y      = ctf.sensor.location(2,index);
          EEG.chanlocs(index).Z      = ctf.sensor.location(3,index);
      end
  catch,
      for index=1:length(ctf.sensor.info),
          EEG.chanlocs(index).labels = ctf.sensor.info(index).label;
          if ~isempty(ctf.sensor.info(index).location)
              EEG.chanlocs(index).X  = ctf.sensor.info(index).location(1);
              EEG.chanlocs(index).Y  = ctf.sensor.info(index).location(2);
              EEG.chanlocs(index).Z  = ctf.sensor.info(index).location(3);
          end;
      end
  end;
  EEG.chanlocs = convertlocs(EEG.chanlocs, 'cart2all');
  
  % now clear the workspace of the input data
  % -----------------------------------------
  clear data
  clear DSpath DSfile DSext i

  % import event information
  % ------------------------
  %try
      eventarray  = [];
      allfields   = {};
      timefields  = {};
      otherfields = {};
      
      eventstruct = ctf_read_markerfile(ctf);
      if isfield(eventstruct, 'markers')
          eventstruct = eventstruct.markers;
      end;
      
      if isfield(eventstruct,'trial_times'),
              
              for index = 1:length(eventstruct)
                  trialtimes = eventstruct(index).trial_times(:,2)*EEG.srate;
                  if EEG.trials > 1
                      trialtimes = trialtimes+(eventstruct(index).trial_times(:,1)-1)*EEG.pnts - ctf.setup.start_msec/1000*EEG.srate;
                  end;
                  
                  for numev = 1:length(trialtimes)
                      EEG.event(end+1).latency = trialtimes(numev);
                      EEG.event(end).type      = eventstruct(index).marker_names;
                      if EEG.trials > 1
                          EEG.event(end).epoch = floor(trialtimes(numev)/EEG.pnts)+1;
                      end;
                  end;
              end;
              
              EEG = eeg_checkset(EEG, 'eventconsistency');
      end
      

  %catch
  %    disp(lasterr);
  %    disp('error (see above) while importing events: events not imported');
  %end;
  
  % command
  com = sprintf('EEG = pop_ctf_read(''%s'', %s)', ctf.folder, vararg2str(options));
    
% format folder
% -------------
function [folder,parentfolder] = formatfolder( orifolder )
    delims = [ find(orifolder == '/') find(orifolder == '\') ];
    if ~isempty(delims)
        if delims(end) == length(orifolder)
            folder = orifolder(delims(end-1)+1:end-1);
            parentfolder = orifolder(1:delims(end-1)-1);
        else
            folder = orifolder(delims(end)+1:end-1);
            parentfolder = orifolder(1:delims(end)-1);
        end;
    else
        folder = orifolder;
        parentfolder = '.';
    end;      
