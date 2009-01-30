function class = ctf_read_classfile(ctf,clsFile);

% class = ctf_read_classfile(ctf,clsFile);
% 
% This function reads the ClassFile.cls created by CTF preprocessing
% software.  
%  
% ctf - a struct returned by ctf_read, it contains ctf.folder, which is a
%   .ds path where a ClassFile.cls can be read; eg:
%   ctf = ctf_read
%   ctf.class = ctf_read_classfile(ctf)
%
% clsFile - this can be used to override the default behavior of reading
%   *.ds/ClassFile.cls, it should be a complete path to a .cls file; eg:
%   ctf.class = ctf_read_classfile([],clsFile)
%   ctf.class = ctf_read_classfile([],'/<completepath>/ClassFile.cls')
%
% It returns a class struct, eg:
%
%   class.path: '/data/dnl3/HighN_cuebase/allTrialData/CJ_HiN_alltrials.ds'
% class.number: 4
%   class.data: {[1x1 struct]  [1x1 struct]  [1x1 struct]  [1x1 struct]}
%
% where class.data is a 1xN struct array, eg:
%
% >> class.data(1)
%
% ans = 
%
%    classgroupid: 3
%            name: 'EYEBLINKa'
%         comment: [1x63 char]
%           color: ''
%        editable: 'Yes'
%         classid: 2
%         Ntrials: 67
%          trials: [1x67 double]
%
% To get all the class names,
% [classNames{1:class.number}] = deal(class.data.name);
% classNames = {ctf.class.data(:).name};
%
% To scan all the class names, use strmatch, eg:
% strmatch('BAD',{ctf.class.data(:).name});
%
% Note that the class.data(n).trials contains numbers beginning with 1,
% whereas the CTF datafile has values beginning with 0 (this is based on
% differences in array indexing between c/c++ and matlab).
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

% Copyright (C) 2005  Darren L. Weber
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

% Created: 12/2005, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~exist('clsFile','var'), clsFile = ''; end
if isempty(clsFile), clsFile = ''; end
if clsFile,
  fid = fopen(clsFile,'r');
else
  fid = -1;
end

if fid < 0,
  % only use the ctf.folder if the above has failed
  if isfield(ctf,'folder'),
    if isempty(ctf.folder),
      error('ctf.folder is empty and no clsFile specified');
    else
      clsFile = fullfile(ctf.folder,'ClassFile.cls');
      fid = fopen(clsFile,'r');
    end
  end
end

if fid < 0,
  error('failed to open any ClassFile.cls');
end


% Read the dataset path from the classfile
class.path = [];
while isempty(class.path),
  curLine = fgetl(fid);
  if strcmp(curLine,'PATH OF DATASET:'),
    curLine = fgetl(fid);
    class.path = curLine;
  end
end
frewind(fid);

% Find the number of classes to read
class.number = [];
while isempty(class.number),
  curLine = fgetl(fid);
  if strcmp(curLine,'NUMBER OF CLASSES:'),
    curLine = fgetl(fid);
    class.number = str2num(curLine);
  end
end
frewind(fid)

if isempty(class.number),
  warning('ClassFile.cls is empty of classes');
  fclose(fid);
  return
end


% read the classes, assuming this format:
%
% CLASSGROUPID:
% 3
% NAME:
% EYEBLINK
% COMMENT:
% Ampl(MEG)=0T; Deriv(MEG)=0T/s; Ampl(EEG)=80uV; Deriv(EEG)=0V/s.
% COLOR:
% 
% EDITABLE:
% Yes
% CLASSID:
% 2
% NUMBER OF TRIALS:
% 3
% LIST OF TRIALS:
% TRIAL NUMBER
%                 +527
%                 +612
%                 +624

% initialize a class.data array of structs
class.data(1:class.number) = struct(...
    'classgroupid', '',...
            'name', '',...
         'comment', '',...
           'color', '',...
        'editable', '',...
         'classid', [],...
         'Ntrials', [],...
          'trials', []);

classCurrent = 0;
while classCurrent < class.number,
  
  curLine = fgetl(fid);
  if strcmp(curLine,'CLASSGROUPID:');
    % update current class no. here only
    classCurrent = classCurrent + 1;
    curLine = fgetl(fid);
    curLine = str2num(curLine);
    class.data(classCurrent).classgroupid = curLine;
    curLine = fgetl(fid);
  end
  
  if strcmp(curLine,'NAME:');
    curLine = fgetl(fid);
    class.data(classCurrent).name = curLine;
    curLine = fgetl(fid);
  end
  
  if strcmp(curLine,'COMMENT:');
    curLine = fgetl(fid);
    class.data(classCurrent).comment = curLine;
    curLine = fgetl(fid);
  end

  if strcmp(curLine,'COLOR:');
    curLine = fgetl(fid);
    class.data(classCurrent).color = curLine;
    curLine = fgetl(fid);
  end

  if strcmp(curLine,'EDITABLE:');
    curLine = fgetl(fid);
    class.data(classCurrent).editable = curLine;
    curLine = fgetl(fid);
  end
  
  if strcmp(curLine,'CLASSID:');
    curLine = fgetl(fid);
    curLine = str2num(curLine);
    class.data(classCurrent).classid = curLine;
    if curLine ~= classCurrent,
      error('curLine ~= classCurrent');
    end
    curLine = fgetl(fid);
  end
  
  if strcmp(curLine,'NUMBER OF TRIALS:');
    curLine = fgetl(fid);
    Ntrials = str2num(curLine);
    class.data(classCurrent).Ntrials = Ntrials;
    
    % read past "LIST OF TRIALS:"
    curLine = fgetl(fid);
    if strcmp(curLine,'LIST OF TRIALS:'),
      % that's good, we expect it
    else
      warning('Did not find ''LIST OF TRIALS:'' where expected')
    end
    
    % read past "TRIAL NUMBER"
    curLine = fgetl(fid);
    if strcmp(curLine,'TRIAL NUMBER'),
      % that's good, we expect it
    else
      warning('Did not find ''TRIAL NUMBER'' where expected')
    end
    
    class.data(classCurrent).trials = [];
    if class.data(classCurrent).Ntrials > 0,
      class.data(classCurrent).trials = zeros(1,Ntrials);
      for n = 1:Ntrials,
        curLine = fgetl(fid);
        % add 1 to trial numbers for matlab compatibility
        class.data(classCurrent).trials(n) = str2num(curLine) + 1;
      end
    end
    
  end

end

return
