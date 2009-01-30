
% ctf2eeglab - script to convert and save ctf .ds into eeglab .set data
%
% ctf2eeglab
%
% The script uses a ctf struct in the matlab workspace or a GUI prompt to
% load a CTF .ds folder, then converts the ctf data into an EEGLAB EEG
% struct, saving the resulting dataset into an EEEGLAB .set file, located
% in the same path as the ctf .ds folder.  The GUI prompt for the CTF .ds
% folder also provides access to definition of the channels, time and
% trials to load.
%

% $Revision: 1.1 $ $Date: 2009-01-30 03:49:26 $

% Licence:  GNU GPL, no express or implied warranties
% Modified: 01/2004, Darren.Weber_at_radiology.ucsf.edu
%                    - developed in collaboration with Fredrick Carver of
%                    the NIH, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.1 $';
fprintf('\nCTF2EEGLAB [v %s]\n\n',ver(11:15));
clear ver

switch exist('eeglab'),
  case 2,
    fprintf('...found eeglab to be a function that exists\n');
  case 7,
    fprintf('...found eeglab to be a directory that exists\n');
  otherwise,
    error('eeglab does not exist');
end

% only prompt for ctf data if it does not already exist in the matlab
% workspace
if exist('ctf') ~= 1,
  
  fprintf('...no ''ctf'' struct in workspace\n');
  [ctf,FIG] = ctf_read_gui;
  uiwait(FIG); clear FIG;
end
if ~isfield(ctf,'data'),
  fprintf('...no ctf.data field in workspace\n');
  [ctf,FIG] = ctf_read_gui;
  uiwait(FIG); clear FIG;
end
if isempty(ctf.data),
  fprintf('...ctf.data field is empty\n');
  [ctf,FIG] = ctf_read_gui;
  uiwait(FIG); clear FIG;
end

% check if the data is averaged
if ctf.setup.number_trials_averaged > 0,
    warning('this .ds folder is averaged');
end

% check if the data is greater than 500 Mb
data_size = ctf.setup.number_samples * ctf.setup.number_channels * ctf.setup.number_trials;
data_bytes = data_size * 8;
if data_bytes > 5e9, warning('data is greater than 500 Mb'); end
clear data_size data_bytes;

% rearrange the trials into a 3D data matrix
% ctf.data is a 3D matrix with samples X channels X trials, 
% whereas EEGLAB is a 3D matrix with channels X samples X trials
clear data
data = zeros( size(ctf.data,2), size(ctf.data,1), size(ctf.data,3) );
for i = 1:size(ctf.data,3),
  data(:,:,i) = ctf.data(:,:,i)';
end


% Startup EEGLAB
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% import the data into the EEGLAB EEG struct

[DSpath,DSfile,DSext] = fileparts(ctf.folder);

EEG.setname = DSfile;
EEG.filename = [DSfile,'.set'];
EEG.filepath = [DSpath,filesep];
EEG.pnts = ctf.setup.number_samples;
EEG.nbchan = ctf.setup.number_channels;
EEG.trials = ctf.setup.number_trials;
EEG.srate = ctf.setup.sample_rate;
EEG.xmin = ctf.setup.start_sec;
EEG.xmax = ctf.setup.end_sec;
EEG.data = data;
EEG.icawinv = [];
EEG.icasphere = [];
EEG.icaweights = [];
EEG.icaact = [];
EEG.event = [];
EEG.epoch = [];
EEG.comments = ctf.setup.run_description';
EEG.ref = 'common';

for i=1:ctf.setup.number_channels,
    EEG.chanlocs(i).labels = ctf.sensor.label{i};
    EEG.chanlocs(i).X      = ctf.sensor.location(1,i);
    EEG.chanlocs(i).Y      = ctf.sensor.location(2,i);
    EEG.chanlocs(i).Z      = ctf.sensor.location(3,i);
end
EEG.chanlocs = pop_chanedit(EEG.chanlocs, 'convert', 'cart2topo', 'convert', 'cart2sph');

% now clear the workspace of the input data
clear ctf data
clear DSpath DSfile DSext i

if ~exist('saveSet','var'), saveSet = 1; end

if saveSet,
  
  [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,...
  'setname', EEG.setname,...
  'comments', EEG.comments,...
  'overwrite','on',...
  'save', [EEG.filepath,filesep,EEG.filename]);
  
else
  
  [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,...
  'setname', EEG.setname,...
  'comments', EEG.comments,...
  'overwrite','on');
  
end

eeglab redraw;
pack
