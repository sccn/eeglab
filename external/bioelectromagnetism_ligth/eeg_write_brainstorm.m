function eeg_write_brainstorm(fileprefix,data)

% eeg_write_brainstorm - Write EEG data into brainstorm format
% 
% Useage: eeg_write_brainstorm(fileprefix, data)
% 
% file  - the path + filename, eg "c:\subj_data.mat"
%         It is actually a fileprefix, as the _data.mat are
%         appended by this function (they are mandatory)
% 
% data  - a matlab structure with brainstorm data fields.
%         There are 2 essential fields (others are initialised):
% 
%         F       a matrix of voltage values (in Volts) with
%                 electrodes in rows and data points in columns.
%         Time    Time is a row vector with time points (in sec)
%                 for each column of the voltage data.
% 
% This script does not convert the voltage data from uV to Volts
% or the Time from msec to sec.  Also, it assumes that the associated
% Channel struct for this EEG data has the last Channel.Type = 'EEG REF'
% so the ChannelFlag array output from this function contains a -1 at 
% index = size(data.F,1) + 1.
% 
% Comment:  See the brainstorm website for more info, at
%           http://neuroimage.usc.edu/.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:54 $

% Licence:  GNU GPL, no express or implied warranties
% Author:   07/2001, Darren.Weber_at_radiology.ucsf.edu
%           10/2002, Darren.Weber_at_radiology.ucsf.edu
%                    modified warnings and channelflag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('EEG_WRITE_BRAINSTORM...\n'); tic;

if isfield(data,'F'),
    F = data.F;
else
    msg = strcat('data structure must contain an ''F'' field.\n',...
                 'F is a matrix of voltage values (in Volts)\n',...
                 'with electrodes in rows and data points in columns.\n');
    error(sprintf(msg));
end
if isfield(data,'Time'),
    Time = data.Time;
else
    msg = strcat('data structure must contain a ''Time'' field.\n',...
                 'Time is a row vector with time points (in sec)\n',...
                 'for each column of the voltage data.\n');
    error(sprintf(msg));
end
if isfield(data,'ChannelFlag'),
    ChannelFlag = data.ChannelFlag;
else
%     msg = strcat('data structure can contain a ''ChannelFlag'' field.\n',...
%                  'ChannelFlag is a row vector indicating the status\n',...
%                  'of the electrodes (-1 is dead, 0 is ignored, 1 is good).\n',...
%                  'Creating this vector = ones(1,size(data.F,1)).\n');
%     fprintf(msg);

    ChannelFlag = [ones(1,size(F,1)), -1];
    fprintf('...assume data Nchannels + 1 is ref, so ChannelFlag(Nchannels + 1) = -1\n');
end
if isfield(data,'NoiseCov'),
    NoiseCov = data.NoiseCov;
else
%     msg = strcat('data structure can contain a ''NoiseCov'' field.\n',...
%                  'This is a square matrix the size of the electrodes\n',...
%                  'Creating this vector = ones(size(data.F,1)).\n');
%     fprintf(msg);
    NoiseCov = [];
end
if isfield(data,'SourceCov'),
    SourceCov = data.SourceCov;
else
%     msg = strcat('data structure can contain a ''SourceCov'' field.\n',...
%                  'This is a square matrix the size of ImageGridLoc(?)\n',...
%                  'Creating this matrix = [].\n');
%     fprintf(msg);
    SourceCov = [];
end
if isfield(data,'Project'),
    Projector = data.Project;
else
%     msg = strcat('data structure can contain a ''Project'' field.\n',...
%                  'Not sure what it is, but it is created for the ascii\n',...
%                  'tutorial data.  Creating this matrix = [].\n');
%     fprintf(msg);
    Project = [];
end
if isfield(data,'Projector'),
    Projector = data.Projector;
else
%     msg = strcat('data structure can contain a ''Projector'' field.\n',...
%                  'This is a matrix of size length of electrodes by rank,\n',...
%                  'not necessarily orthogonal due to Good-Channel selections.\n',...
%                  'Orthogonalized as ''U'', then the data are to be projected away\n',...
%                  'from the projector as F'' = F - U*(U^t * F).\n',...
%                  'Creating this matrix = [].\n');
%     fprintf(msg);
    Projector = [];
end
if isfield(data,'Comment'),
    Comment = data.Comment;
else
%     msg = strcat('data structure can contain a ''Comment'' field.\n',...
%                  'This is a character string describing the data.\n',...
%                  'Creating this field = file.\n');
%     fprintf(msg);
    Comment = fileprefix;
end
if isfield(data,'Device'),
    Device = data.Device;
else
%     msg = strcat('data structure can contain a ''Device'' field.\n',...
%                  'This is a character string describing the device\n',...
%                  'used to acquire the data (eg, SYNAMPS)\n',...
%                  'Creating an empty field.\n');
%     fprintf(msg);
    Device = 'Neuromag_Planar_122'; % the only supported device!
end

[path,file,ext] = fileparts(fileprefix);
ext = '_data.mat';

file = fullfile(path,[file ext]);

save(file,'F','Time','ChannelFlag','NoiseCov','SourceCov',...
          'Project','Projector','Comment','Device');

fprintf('...wrote brainstorm data to:\n\t%s\n', file);

t = toc; fprintf('...done (%6.2f sec).\n\n',t);

return
