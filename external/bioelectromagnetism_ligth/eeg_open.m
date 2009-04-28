function [p] = eeg_open(p,parent)

% eeg_open - function to handle various eeg_load commands
% 
% Usage: [p] = eeg_open(p,[parentgui])
% 
% p is a parameter structure. See eeg_toolbox_defaults for more 
% information on this parameter structure.
% 
% In this function, p must contain the fields:
% 
% p.volt.path - the directory location of the file to load
% p.volt.file - the name of the file to load
% p.volt.type - the file format string, one of:
% 
% 'ASCII'
% 'EMSE'
% 'Scan4x'
% 'Scan3x'
% 'Matlab'
%
% These are the only file types currently supported. 
% See functions eeg_load* for details.
% 
% The most important return value is the ERP data in 
% p.volt.data.  If the file format is scan4x, various 
% ERP parameters are returned also.
% 
% If the ASCII or Matlab type is given, the routine will try 
% to load an associated variance file.  This file must be 
% located in the same path, with the same file name as the 
% voltage file, but the file extension should be '.var'
% 
% See also: EEG_LOAD, EEG_LOAD_ASCII, EMSE_READ_AVG, 
%           EEG_LOAD_SCAN4_AVG, EEG_LOAD_SCAN3_AVG
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:52 $

% Licence:  GNU GPL, no express or implied warranties
% History:  02/2002, Darren.Weber_at_radiology.ucsf.edu
%           04/2002, Darren.Weber_at_radiology.ucsf.edu
%                    added variance handling
%           08/2002, Darren.Weber_at_radiology.ucsf.edu
%                    added EMSE avg handling
%                    added interpolation of zero point option
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('p','var'),[p] = eeg_toolbox_defaults; end

eegversion = '$Revision: 1.1 $';
fprintf('EEG_OPEN [v %s]\n',eegversion(11:15)); tic;


[path,name,ext] = fileparts(strcat(p.volt.path, filesep, p.volt.file));
file = fullfile(path,[name ext]);

if ~isequal(exist(file),2),
  lookfile = which(file);
  if isempty(lookfile),
    msg = sprintf('...cannot locate %s\n', file);
    error(msg);
  else
    file = lookfile;
  end
end

type = lower(p.volt.type)

switch type,
  
  case 'ascii',
    
    [ p.volt.data, p.volt.var ] = eeg_load_ascii(file);
    
  case 'emse',
    
    avg = emse_read_avg(file);
    
    % load variance data?
    if isfield(avg,'volt'),
      if isempty(avg.volt),
        error('failed to read file.');
      end
    else
      error('failed to read file.');
    end
    
    p.volt.data          = avg.volt;
    %p.volt.var           = avg.variance;
    %p.volt.channelNames  = avg.chan_names;
    p.volt.channels      = avg.channels;
    p.volt.points        = avg.pnts;
    p.volt.sampleMsec    = avg.rate;
    p.volt.epochStart    = avg.xmin;
    %p.volt.epochEnd      = avg.xmax;
    p.volt.sampleHz      = 1000 / avg.rate;
    %p.volt.sweeps        = avg.nsweeps;
    
    
  case 'scan4x_cnt',
    
    p.cnt = eeg_load_scan4_cnt(file);
    
  case 'scan4x',
    
    avg = eeg_load_scan4_avg(file);
    
    if ~isfield(avg,'data'),
      msg = sprintf('...failed to load scan4.x avg file:\n... %s\n',file);
      error(msg);
    end
    if isequal([avg.header.domain],1),
      msg = sprintf('...cannot open frequency domain file:\n... %s\n',file);
      error(msg);
    end
    
    p.volt.points = avg.header.pnts;
    p.volt.sampleHz = avg.header.rate;
    p.volt.sampleMsec = 1000/ avg.header.rate;
    p.volt.channels = avg.header.nchannels;
    p.volt.epochStart = avg.header.xmin * 1000; % convert to msec
    p.volt.epochEnd   = avg.header.xmax * 1000; % convert to msec
    p.volt.sweeps = avg.header.acceptcnt;
    
    ampData = [avg.data.samples]; % elect in col, samples in row
    baseline = repmat([avg.electloc.baseline],p.volt.points,1);
    calibration = repmat([avg.electloc.calib],p.volt.points,1);
    n = repmat([avg.electloc.n],p.volt.points,1);
    % Convert to uV
    p.volt.data = ( ampData - baseline ) .* calibration ./ n;
    p.volt.var = [avg.variance.samples];
    
    
  case 'scan3x',
    
    avg = eeg_load_scan3_avg(file);
    
    if isfield(avg,'signal'),
      p.volt.data          = avg.signal;
      p.volt.var           = avg.variance;
      p.volt.channelNames  = avg.chan_names;
      p.volt.points        = avg.pnts;
      p.volt.sampleHz      = avg.rate;
      p.volt.epochStart    = avg.xmin;
      p.volt.epochEnd      = avg.xmax;
      p.volt.sampleMsec    = avg.rate / 1000;
      p.volt.sweeps        = avg.nsweeps;
    else
      msg = sprintf('...failed to load scan3.x datafile: %s\n',file);
      error(msg);
    end
    
  case 'matlab',
    
    p.volt.data = eeg_load(file);
    % Attempt to load associated variance file
    varfile = fullfile(path,strcat(name,'.var'));
    if isequal(exist(varfile),2),
      p.volt.var = eeg_load(varfile);
    else
      lookfile = which(varfile);
      if isempty(lookfile),
        fprintf('...cannot locate Matlab variance file:\n... %s\n', varfile);
      else
        varfile = lookfile;
        p.volt.var = eeg_load(varfile);
      end
    end
    
  otherwise,
    msg = sprintf('\nPlease specify voltage data type: ASCII | EMSE | Scan4x | Scan3x | Matlab?\n\n');
    error(msg);
end

% Try to arrange electrodes in columns (assuming more sample points than electrodes)
s = size(p.volt.data);
if s(1) < s(2),
  fprintf('...rotating voltage data from %s : ',[name ext]);
  p.volt.data = p.volt.data';
  s = size(p.volt.data);
  fprintf('%d rows, %d cols\n', s(1), s(2));
end
p.volt.channels = s(2);
p.volt.points = s(1);

s = size(p.volt.var);
if s(1) < s(2),
  if isequal(p.volt.type,'ASCII'), file = varfile; end
  if isequal(p.volt.type,'Matlab'), file = varfile; end
  fprintf('...rotating variance data from:\n... %s : ',varfile);
  p.volt.var = p.volt.var';
  s = size(p.volt.var);
  fprintf('%d rows, %d cols\n', s(1), s(2));
end


% Verify that essential ERP parameters are set
if isempty(p.volt.sampleHz) | isempty(p.volt.epochStart) | isempty(p.volt.epochEnd),
  
  if exist('parent','var'),
    if ~isempty(parent),
      %help = helpdlg('Please specify ERP parameters as follows...','EEG OPEN HELP');
      %movegui(help,'center'); waitfor(help);
      
      data = get(parent,'UserData');
      data.p = p;
      set(parent,'UserData',data);
      
      tmpgui = gui_eeg_ascii_parameters(parent);
      
      data = get(parent,'UserData');
     [p] = data.p; 
      clear data parent tmpgui;
    end
  else
    fprintf('...please specify ERP parameters in p structure.\n');
  end
end


% --- Setup data structures for timing array

% Interpolate the Zero value
if p.volt.interpZero,
  
  if p.volt.epochStart & p.volt.epochEnd & p.volt.sampleMsec,
    
    if mod(p.volt.epochStart,fix(p.volt.epochStart)) < .1,
      % round the epochStart value, if its remainder is within .1 msec
      p.volt.epochStart = fix(p.volt.epochStart);
    end
    
    p.volt.timeArray = [p.volt.epochStart:p.volt.sampleMsec:p.volt.epochEnd]';
    
    timeNonZero = find(p.volt.timeArray);
    timeZero = find(p.volt.timeArray == 0);
    
    volt = p.volt.data;
    
    InterpZero = interp1( p.volt.timeArray(timeNonZero), volt, 0, 'cubic' );
    volt = [volt(1:timeZero-1,:); InterpZero; volt(timeZero:end,:)];
    p.volt.data = volt; 
    
    clear InterpZero timeNonZero timeZero volt;
  end
else
  p.volt.timeArray = [1:p.volt.points]';
  if and(p.volt.epochStart,p.volt.sampleMsec),
    t = p.volt.epochStart;
    for i=1:p.volt.points,
      p.volt.timeArray(i,1) = t;
      t = t + p.volt.sampleMsec;
    end
  end   
end

p.volt.points = size(p.volt.data,1);

t = toc; fprintf('...done (%6.2f sec).\n\n',t);

return
