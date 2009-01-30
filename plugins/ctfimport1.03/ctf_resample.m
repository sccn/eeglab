function ctf = ctf_resample(folder,newFreq)
% ctf_resample - Resamples a CTF .ds folder to a different sampling rate.
%
% ctf = ctf_resample(folder,newFreq)
%
% The function modifies two files in the .ds folder:
%
% *.meg4 - the raw data
% *.res4 - the data resource file (essentially a header). 
%
% This function calls the matlab resample function.  There is also a
% -resample option in the ctf software command, newDs, but it does not
% specify the sample rate in Hz.
%
% INPUTS: folder (optional - brings up a file selector if not specified):
%         This is the .ds folder containing both the .meg4 and the .res4
%         files
% NEWFREQ: Frequency (in Hz) that you want to pretend the data were sampled
%          at. This can be lower or higher than the actual sampling
%          frequency.
% RETURNS: 'ctf': the ctf data structure (without the raw data): contains the new
%           parameters for sample_rate etc 
%

% $Revision: 1.1 $ $Date: 2009-01-30 03:49:27 $
%
% Copyright (C) 2003  Alex Wade wade@ski.org
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

% Modified: 05/2003, Alex Wade wade@ski.org
%                    - Wrote it
%           Not tested with pretrigger points or averaged data sets.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%--------------------------------------------------------------
ver = '$Revision: 1.1 $';

fprintf('\nCTF_RESAMPLE [v %s]\n',ver(11:15)); tic;

% Populate the 'ctf' data structure
if exist('folder','var')
    if(~isempty(folder))
        ctf = ctf_read_meg4(folder,[],[],[0 0]); % Read in no actual data
    else
        ctf = ctf_read_meg4([],[],[],[0 0]); % 
    end
else
    ctf = ctf_read_meg4([],[],[],[0 0]); % 
end






if ~exist('COEFS','var'),
    COEFS = false;
end

if ~exist('folder','var'),
  if ~exist('ctf','var'),
    ctf = ctf_folder;
  else
    ctf = ctf_folder([],ctf);
  end
else
  if ~exist('ctf','var'),
    ctf = ctf_folder(folder);
  else
    ctf = ctf_folder(folder,ctf);
  end
end

if ~isfield(ctf,'setup'),
  ctf = ctf_read_res4(ctf.folder,1,COEFS);
end





% RESAMPLING GOES HERE...
% A typical ctf data structure looks like...
%      folder: [1x97 char]
%      res4: [1x1 struct]
%      setup: [1x1 struct]
%      sensor: [1x1 struct]
%      meg4: [1x1 struct]
%      data: [294000x275 double]

% The fields to change are:
% ctf.setup.sample_rate
% ctf.setup.number_samples
% ctf.setup.sample_sec
% ctf.setup.sample_msec
% ctf.setup.pretrigger_samples (if there were any)
% ctf.setup.run_description (to append the fact that we've resampled the
% data)
% 
% 

oldFreq = ctf.setup.sample_rate;

nChannels = ctf.setup.number_channels;
nPoints = ctf.setup.sample_rate * ctf.setup.duration_trial;
nTrials = ctf.setup.number_trials;
newNPoints = round(nPoints * newFreq / oldFreq); % This should match 'resample's idea of the new number of points

ctf.setup.sample_rate = newFreq;
ctf.setup.number_samples = newNPoints;
ctf.setup.sample_sec = (ctf.setup.sample_sec * oldFreq / newFreq);
ctf.setup.sample_msec = ctf.setup.sample_sec * 1000;
ctf.setup.pretrigger_samples = round(ctf.setup.pretrigger_samples * newFreq / oldFreq);

% When we write out the meg4 file, we will just skip the first 8 bytes (the
% header 'MEG41CP'+13? and then start fwriting (signed) 4-byte integers
% ('int32')

fidOld = fopen(ctf.meg4.file,'rb+','s'); % 'ieee-be.l64' or 's' - IEEE floating point with big-endian byte
fidNew = fopen([ctf.meg4.file,'.resamp'],'wb','s'); 
% ordering and 64 bit long data type.

if (fidOld<0)
    error('Could not open meg4 file for reading - does it exist?');
end
if (fidNew<0)
    fclose(fidOld);
    error('Could not open new meg4 file for writing - check file permissions');
end

disp('Copying header');
header=fread(fidOld,8,'char'); % Start writing at the 9th byte (first 8 are the header)
cnt=fwrite(fidNew,header,'char');
if (cnt~=8)
    fclose(fidOld);
    fclose(fidNew);
    error('Failed to write to the new meg4 file.');
end
% Loop over trials and channels
disp('Resampling...');
h=waitbar(0,'Resampling');

% Resample will work on the columns of a data matrix. But it runs out of
% memory if we try to process the whole thing at once. So we will do it
% channel by channel at the same time that we write it out...

for thisTrial = 1:nTrials
    for thisChannel = 1:nChannels
        % Do resampling on raw data from disk - not stuff stored in
        % ctf.data. Avoids any messing around with channel gains...
        
        origData = fread(fidOld,nPoints,'int32');
        % Using the default matlab resamp params. 
        newDataVector = round((resample(double(origData),newFreq,oldFreq)));
        
        cnt = fwrite(fidNew,newDataVector,'int32');
        
        if (cnt ~= newNPoints) % Check for write success
            disp(cnt);
            disp(newNPoints);
            fclose(fidNew);
            error('Failed to write all the data - check to see if the disk is full and whether you have the correct permissions');
        end
        waitbar(thisChannel/nChannels);    
    end
end

close(h); % Close the progress bar

% Close the files
fclose(fidOld);
fclose(fidNew);
fprintf('\nDone meg4\n');

% Finished writing the raw data....

% Now create a new .res4 file. 
fidOld = fopen(ctf.res4.file,'rb','ieee-be.l64');
fidNew = fopen([ctf.res4.file,'.resamp'],'wb','ieee-be.l64');

if fidOld < 0, error('cannot open original.res4 file'); end
if fidNew < 0
    fclose(fidOld);
    error('Cannot open temp.res4 file'); 
end

% Copy the first set of data from the file
buff = fread(fidOld,1288,'uint8');
cnt = fwrite(fidNew,buff,'uint8');
if (cnt ~= length(buff))
    fclose(fidOld);
    fclose(fidNew);
    error('Failed to write to the new res4 file.');
end

%---NUMBER OF SAMPLES
dummy = fread(fidOld,1,'int32');
cnt = fwrite(fidNew,ctf.setup.number_samples,'int32');

% Only check writing once
if (cnt~=1)
    fclose(fidOld);
    fclose(fidNew);
    error('Failed to write to the new res4 file.');
end

%---SAMPLE RATE
% Copy next 4 bytes
buff = fread(fidOld,4,'uint8');
cnt = fwrite(fidNew,buff,'uint8');
if (cnt ~= length(buff))
    fclose(fidOld);
    fclose(fidNew);
    error('Failed to write to the new res4 file.');
end
dummy = fread(fidOld,1,'double');
cnt = fwrite(fidNew,ctf.setup.sample_rate,'double');
if (cnt ~= 1)
    fclose(fidOld);
    fclose(fidNew);
    error('Failed to write to the new res4 file.');
end
% copy next 16 bytes
buff = fread(fidOld,16,'uint8');
cnt = fwrite(fidNew,buff,'uint8');
if (cnt ~= length(buff))
    fclose(fidOld);
    fclose(fidNew);
    error('Failed to write to the new res4 file.');
end

%---PRETRIGGER POINTS
cnt = fwrite(fidNew,ctf.setup.pretrigger_samples,'int32');
if (cnt ~= 1)
    fclose(fidOld);
    fclose(fidNew);
    error('Failed to write to the new res4 file.');
end

dummy = fread(fidOld,1,'int32');

buff = fread(fidOld,inf,'uint8');
cnt = fwrite(fidNew,buff,'uint8');
if (cnt ~= length(buff))
    fclose(fidOld);
    fclose(fidNew);
    error('Failed to write to the new res4 file.');
end

fclose(fidOld);
fclose(fidNew);


return
