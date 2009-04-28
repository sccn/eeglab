function [p] = eeg_toolbox_defaults(command,p);

% eeg_toolbox_defaults - Create, read, write eeg_toolbox defaults
%
% Useage: [p] = eeg_toolbox_defaults(command,[p])
%
% command  =  'create'
%             'read'
%             'save'   | 'write'
%             'saveas' | 'write_other'
% p        =  structure returned by create|read command
%             structure to write for write command
% 
% All read and write commands are to a matlab .mat file.  The
% 'saveas' or 'write_other' command will write out a data specific archive
% that can be accessed from the recent files list.  It will be located in
% the folder where the eeg data is opened.
% 
% Examples:
%[p] = eeg_toolbox_defaults; % create new defaults
%[p] = eeg_toolbox_defaults('read'); % read saved defaults
% eeg_toolbox_defaults('save',p);  % write current p as default
% eeg_toolbox_defaults('saveas',p); % write current p data
%
% The save or write command will write to the eeg_toolbox
% installation folder.  The read command will first 
% try to find the paramater file in the eeg_toolbox
% installation and otherwise recreates the defaults.
%
% Notes:    Handles parameter structure for the 
%           eeg_toolbox routines.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:53 $

% Licence:  GNU GPL, no express or implied warranties
% Created:  01/2002 Darren.Weber_at_radiology.ucsf.edu
% Modified: 02/2002 Darren.Weber_at_radiology.ucsf.edu
%           - adapted read/write output from .txt format 
%             to .mat format so it will contain data structures.
%             Hence, it is no longer editable text.
%           08/2002 Darren.Weber_at_radiology.ucsf.edu
%                   added MRI defaults
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eegversion = '$Revision: 1.1 $';
fprintf('\nEEG_TOOLBOX_DEFAULTS [v %s]\n',eegversion(11:15));

if version('-release') < 10,
    msg = printf('matlab release < version 10, eeg_toolbox GUIs may fail.\n');
    warning(msg);
end

if ~exist('command','var'), command = 'create'; end

% try to locate the installation path
eegPath = eeg_toolbox_path;


% try to locate the mri_toolbox installation path
mriPath = fileparts(which('avw_view'));
if isempty(mriPath),
    msg = sprintf('Cannot find mri_toolbox on the matlab path.\nPlease install and use the addpath command.\n\n');
    warning(msg);
    mriPath = eegPath;
else
    mriPath = strcat(mriPath,filesep);
end

% remove any saved preferences
if ispref('eeg_toolbox'),
    rmpref('eeg_toolbox');
end


switch command
case 'create',
    
    fprintf('...creating eeg_toolbox defaults\n\n');
    
   [p] = [];
    
    % default location for parameter file
    p.file = strcat(eegPath,'eeg_example_data',filesep,'eeg_toolbox_defaults.mat');
    
    % parameters for 'elec_load' routine, see that for more information
    p.elec.path     = strcat(eegPath,'eeg_example_data',filesep);
    p.elec.file     = 'elec_124_cart.txt';
    p.elec.n        = 129;         % how many electrodes to load
    p.elec.type     = 'Cartesian'; % see elec_load for avail. options
    p.elec.shape    = 'sphere';    % original = '', 'sphere', or 'ellipse'
    p.elec.data     = [];
    p.elec.plot     = 0; % plot electrode positions
    p.elec.plotSurf = 1; % plot electrode surface (from convhull)
    
    % parameters for 'mesh_open' & 'mesh_plot', see these for more information
    p.mesh.path         = strcat(eegPath,'eeg_example_data',filesep);
    p.mesh.file         = 'BrainStorm_subjecttess.mat';
    p.mesh.type         = 'BrainStorm'; % see mesh_open.m for options
    p.mesh.data         = []; % init for mesh_open.m
    p.mesh.plot.open    = 0; % boolean arg for mesh_open.m
    p.mesh.plot.overlay = 0; % Overlay arg in mesh_plot.m
    p.mesh.plotSurf     = 0; % plot scalp mesh surface (eeg_contours_engine)
    p.mesh.reload       = 1; % boolean, 1 = reload mesh data
    p.mesh.samplePoint  = [];
    p.mesh.sampleTime   = []; % Time of current sample point (msec)
    p.mesh.sampleMsec   = []; % Sample rate (msec)
    p.mesh.sampleHz     = []; % Sample rate (Hz)
    
    % parameters for 'eeg_open' routine, see that for more information
    p.volt.path                 = strcat(eegPath,'eeg_example_data',filesep);
    p.volt.file                 = 'eeg_124ch.dat';
    p.volt.type                 = 'ascii'; % see eeg_open.m for valid options
    
    p.volt.data                 = []; % potential (uV)
    p.volt.var                  = []; % variance for ERP (uV^2)
    p.volt.timeArray            = []; % The timing of each sample point (msec)
    p.volt.points               = []; % Number of sample points in epoch
    p.volt.sampleHz             = 400; % Sample rate (Hz)
    p.volt.sampleMsec           = 2.5; % Sample rate (msec)
    p.volt.sampleTime           = 400; % Time of current sample point (msec)
    p.volt.samplePoint          = 250; % single or start time point (row of p.volt.data)
    p.volt.channels             = []; % Number of channels/electrodes
    p.volt.epochStart           = -200; % Time of epoch start (msec)
    p.volt.epochEnd             = 1500; % Time of epoch end (msec)
    p.volt.sweeps               = []; % Number of accepted sweeps
    p.volt.peaks                = []; % Data array of ERP peak values
    p.volt.interpZero           = 1;  % Interpolate time zero ERP potential,
    % mainly for ascii versions of NeuroScan data
    
    p.cnt.path                  = strcat(eegPath,'eeg_example_data',filesep);
    p.cnt.file                  = 'scan41_short.cnt';
    p.cnt.type                  = 'scan4x_cnt';
    
    p.eeg.path                  = strcat(eegPath,'eeg_example_data',filesep);
    p.eeg.file                  = 'scan41_short.eeg';
    p.eeg.type                  = 'scan4x_eeg';
    
    
    
    % parameters for 'erf_open' routine, see that for more information
    p.erf.path                 = strcat(eegPath,'eeg_example_data',filesep);
    p.erf.file                 = 'eeg_124ch.dat';
    p.erf.type                 = 'ascii'; % see eeg_open.m for valid options
    
    p.erf.data                 = [];   % potential (Tesla)
    p.erf.var                  = [];   % variance for ERF (T^2)
    p.erf.timeArray            = [];   % The timing of each sample point (msec)
    p.erf.points               = [];   % Number of sample points in epoch
    p.erf.sampleHz             = 400;  % Sample rate (Hz)
    p.erf.sampleMsec           = 2.5;  % Sample rate (msec)
    p.erf.sampleTime           = 400;  % Time of current sample point (msec)
    p.erf.samplePoint          = 250;  % single or start time point (row of p.erf.data)
    p.erf.channels             = [];   % Number of channels/electrodes
    p.erf.epochStart           = -200; % Time of epoch start (msec)
    p.erf.epochEnd             = 1500; % Time of epoch end (msec)
    p.erf.sweeps               = [];   % Number of accepted sweeps
    p.erf.peaks                = [];   % Data array of ERP peak values
    p.erf.interpZero           = 0;    % Interpolate time zero ERF potential
    
    
                                       
    % parameters for 'eeg_contours_engine', see that for more information
    p.rangeMethod                 = 'minmaxabs' ; % minmax? = 'abs','one','all','giv'
    p.minimumIntensity            =  -10    ;	
    p.maximumIntensity            =   10    ;
    p.interpMethod                = 'cubic' ;   % matlab griddata interpolation method
    
    p.timeMethod                  =    1    ;   % 1=single time point, 2=range of time points
    p.endTime                     =    2    ;   % finish time point (row of p.volt.data)
    
    p.contour.stepMethod          =  0; % 0 is step size, 1 is number of steps
    p.contour.Nsteps              = 10;
    p.contour.stepSize            =  1; % 1 uV or similar units
    p.contour.raw2D               =  0; % plot 2D contours
    p.contour.plot2D              =  0; % plot projected 2D contours
    p.contour.plot3D              =  0; % plot 3D surface contours
    
    % Now defunct (02/2003), as topo maps now use triangulations, rather than meshgrid
    % keeping for now, just in case they prove useful at some later date
    p.grid.method                  =    1    ;   % Use grid size method
    p.grid.size                    =  100    ;
    p.grid.sizeMin                 =   10    ;
    p.grid.sizeMax                 =  210    ;
    p.grid.res                     =    0.1  ;   % 1 mm
    p.grid.resMin                  =    0.01 ;
    p.grid.resMax                  =    2.01 ;
    
    % Miscellaneous
    p.clickTimePoint              = 0;   % select a time point interactively
    p.saveGraphics                = 'No Save Plots'; % Do not save graphic files by default
    p.topoView                    = '2D'; % Default topography view
    p                             = eeg_colormap(p); % Return defaults for 'eeg_colormap'
    p.hold                        = 0;   % default for GUI hold checkboxes
    
case 'read'
    
    % Look for default file 
    [path,name,ext] = fileparts(strcat(eegPath,'eeg_toolbox_defaults.mat'));
    p.file = fullfile(path,[name ext]);
    file = exist(p.file);
    if file == 2,
        fprintf('...reading eeg_toolbox defaults from:\n%s\n\n',p.file);
        load(p.file);
    else
        fprintf('...cannot locate eeg_toolbox defaults - recreating defaults.\n\n');
       [p] = eeg_toolbox_defaults;
    end
    
    % verify that path to default files exists;
    % if not, create defaults again (I hope this will avoid 
    % new installation problems)
    voltPathExist = exist(p.volt.path);
    if ~voltPathExist,
        fprintf('...voltage path does not exist - reinitializing defaults.\n');
       [p] = eeg_toolbox_defaults('create');
        eeg_toolbox_defaults('write',p);
    end
    
case {'write','save'}
    
    if ~exist('p','var'),
        msg = sprintf('p argument is essential: eeg_toolbox_defaults(''write'',p)');
        error(msg);
    end
    
    [path,name,ext] = fileparts(strcat(eegPath,'eeg_toolbox_defaults.mat'));
    p.file = fullfile(path,[name ext]);
    
    save(p.file,'p');
    
    fprintf('...saved eeg_toolbox defaults to:\n%s\n\n',p.file);
    
case {'write_other','saveas'}
    
    if ~exist('p','var'),
        msg = sprintf('p argument is essential: eeg_toolbox_defaults(''write_other'',p)');
        error(msg);
    end
    
    [path,name,ext] = fileparts(strcat(p.volt.path, filesep, p.volt.file));
    ext = '.mat';
    p.file = fullfile(path,[name ext]);
    
    
    [newfile,newpath] = uiputfile(p.file,'Save EEG_TOOLBOX Parameters to File:');
    if (newfile == 0) & (newpath == 0),
        fprintf('...aborting save\n\n');
    else
        [path,name,ext] = fileparts(strcat(newpath, filesep, newfile));
        ext = '.mat';
        p.file = fullfile(path,[name ext]);
        save(p.file,'p');
        fprintf('...saved eeg_toolbox parameters to:\n%s\n\n',p.file);
        recentfiles = eeg_toolbox_recent(p.file);
    end
    
otherwise
    error('...invalid command\n\n');
end

return
