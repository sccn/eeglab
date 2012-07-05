% create_study_fv() - creates a STUDY structure using ALL sets stored in a folder specified by the parameter "pathin".
%
% Usage: 
%             >> [STUDY, ALLEEG] = create_study_fv(pathin,pathout,studyname,badcomps);  
% Inputs:
%   pathin     - [string] path to the folder containing datasets
%   pathout    - [string] path where STUDY is going to be saved
%   studyname  - [string] name used to save STUDY
%   badcomps   - [0 1] badcomps=1 -> EEG.badcomps is created for each
%                single dataset as an empty field; badcomps=0 -> EEG.badcomps is not created 
%
% Outputs:
%   STUDY      - a new STUDY set containing ALL datasets stored in the folder 
%   ALLEEG     - a vector of EEG datasets included in the STUDY structure 
%
%  See also:  pop_createstudy()
%
% Authors: Filipa Campos Viola, 24/01/2008, MRC-IHR, Southampton, UK
% (f.viola@soton.ac.uk)

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) F. Campos Viola, MRC-IHR
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

% revised by F Campos-Viola - corrmap1.01 (30/01/2009)

% revised by F Campos-Viola - adding option to created EEG.badcomps as an empty field for each dataset (30/01/2009)

function [STUDY,ALLEEG]=create_study_fv(pathin, pathout,studyname,badcomps)

cd(pathin) %points to a path
list=dir('*.set'); %reads all .set files in that path


%[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

STUDY=[];
ALLEEG=[];

n=length(list); %number of datasets that are going to be included in STUDY

if badcomps==1
for j=1:n
    
    EEG = pop_loadset(list(j).name, pathin);
    
    EEG.badcomps=[];
    
    EEG = eeg_checkset(EEG);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG);
    
       
    EEG = pop_saveset(EEG, list(j).name, pathin);

end
end

for i= 1:n

    %creating the study
    [STUDY ALLEEG]= std_editset(STUDY, ALLEEG,'name', studyname, 'commands', {...
        {'index', i,'load', [pathin,list(i).name],'subject',num2str(i)}
    });
end

%saving the STUDY and not changing datasets
[STUDY ALLEEG]= std_editset(STUDY,ALLEEG,'filename',[pathout studyname '.study'], 'updatedat', 'off');

%CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];

[STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG);

