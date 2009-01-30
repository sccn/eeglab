function ctf_write_ascii(OutputFilePrefix,ctf,CHAN,TIME,TRIALS)

% ctf_write_ascii - write ctf.data into ascii file
%
% ctf_write_ascii(OutputFilePrefix,ctf,CHAN,TIME,TRIALS)
%
% OutputFilePrefix - a file prefix is used so that data information can be
% appended to the file name.  The default value is to use the ctf.folder
% filename.
%
% CHAN - see ctf_channel_select for options
% TIME - see ctf_read for options (given in msec)
% TRIALS - select 1 trial to plot (the default is trial = 1)
%
% The output data matrix is time (rows) x channels (columns)
%
%      <>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%      <                                                      >
%      <                    DISCLAIMER:                       >
%      <                                                      >
%      < THIS PROGRAM IS INTENDED FOR RESEARCH PURPOSES ONLY. >
%      < THIS PROGRAM IS IN NO WAY INTENDED FOR CLINICAL OR   >
%      <                    OFFICIAL USE.                     >
%      <                                                      >
%      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>
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

% Modified: 07/2004, Darren.Weber_at_radiology.ucsf.edu
%                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.1 $';
fprintf('\nCTF_WRITE_ASCII [v %s]\n',ver(11:15)); tic;

if ~exist('ctf','var'),
    error('no input ctf data struct!');
end

if ~exist('OutputFilePrefix','var'),
    [ctfPath, ctfFolder, ext] = fileparts(ctf.folder);
    OutputFilePrefix = ctfFolder;
end
if isempty(OutputFilePrefix),
    [ctfPath, ctfFolder, ext] = fileparts(ctf.folder);
    OutputFilePrefix = ctfFolder;
end

if ~exist('CHAN','var'), CHAN = 'all'; end
if ~exist('TIME','var'), TIME = 'all'; end
if ~exist('TRIALS','var'), TRIALS = 1; end

if isempty(CHAN), CHAN = 'all'; end
if isempty(TIME), TIME = 'all'; end
if isempty(TRIALS), TRIALS = 1; end

% This function calls ctf_channel_select
[CHAN,type] = ctf_channel_select(ctf,CHAN);

switch num2str(TIME),
  case 'all',
    TIME = ctf.setup.time_msec;
    TIME_index = 1:ctf.setup.number_samples;
otherwise
    fprintf('...sorry, time restrictions not implemented correctly (03/2004)\n');
    TIME = ctf.setup.time_msec;
    TIME_index = 1:ctf.setup.number_samples;
    
    
    %     % assume the input is a range of times in sec
    %     % check the range
    %     if TIME(1) > ctf.setup.time_msec(1),
    %       fprintf('...setting TIME(1) = ctf.setup.time_msec(1)\n');
    %       TIME(1) = ctf.setup.time_msec(1);
    %     end
    %     if TIME(end) > ctf.setup.time_msec(end),
    %       fprintf('...setting TIME(end) = ctf.setup.time_msec(end)\n');
    %       TIME(end) = ctf.setup.time_msec(end);
    %     end
    %     % now find the nearest indices into the samples matrix
    %     TIME_index = interp1(ctf.setup.time_msec,1:ctf.setup.number_samples,TIME,'nearest');
    %     % now ensure that the TIME array is consistent with ctf.setup.time_sec
    %     TIME = ctf.setup.time_msec(TIME_index);
end
TIME = sort(TIME);



Ntrials = size(ctf.data,3);
switch num2str(TRIALS),
  case 'all',
    TRIALS = 1:ctf.setup.number_trials;
  otherwise
    % assume the input is an array of trials
end
TRIALS = unique(sort(TRIALS));



for trial = TRIALS,
    
    fprintf('...exporting trial %d of %d trials in ctf.data\n',trial,Ntrials);
    
    % now export the data
    
    data = ctf.data(TIME_index,CHAN,trial);
    
    OutputFileName = [OutputFilePrefix,'_trial',num2str(trial),'.txt'];
    
    save(OutputFileName, 'data', '-ascii', '-double', '-tabs')
    
end


t = toc; fprintf('...done (%6.2f sec)\n\n',t);

return
