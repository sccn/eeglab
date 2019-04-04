% std_limoerase() - erase LIMO files.
%
% Usage:
%   >> std_limoerase(filepath, filename, subjects, designname)
%
% Inputs:
%   filepath     - [string] STUDY file path
%   filename     - [string] STUDY file name
%   subjects     - [cell] list of subject name
%   designname   - [string] name or number of the design
%
% Author: Ramon Martinez and Arnaud Delorme, SCCN, UCSD, 2018-
%
% See also: std_limo()

% Copyright (C) Ramon Martinez and Arnaud Delorme
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function std_limoerase(filepath, filename, subjects, designname)

if nargin < 4
    help std_limoerase;
    return;
end

if exist([filepath filesep 'limo_batch_report'],'dir')
    try
        rmdir([filepath filesep 'limo_batch_report'],'s');
    catch
    end
end

% Cleaning level 1 folders
for i = 1:length(subjects)
    tmpfiles = dir([filepath filesep 'LIMO_' filename filesep subjects{i} filesep 'GLM' designname '*']);
    tmpfiles = {tmpfiles.name};
    for j = 1:length(tmpfiles)
        try
            rmdir([filepath filesep 'LIMO_' filename filesep subjects{i} filesep tmpfiles{j}],'s');
        catch
        end
    end
end

% Cleaning level 2 folders
% if isfield(STUDY.design(design_index),'limo') && isfield(STUDY.design(design_index).limo,'groupmodel')
%     for nlimos = 1:length(STUDY.design(design_index).limo)
%         for ngroupmodel = 1:size(STUDY.design(design_index).limo(nlimos).groupmodel,2)
%             path2file = fileparts(STUDY.design(design_index).limo(nlimos).groupmodel(ngroupmodel).filename);
%             try
%                 rmdir(path2file,'s');
%             catch
%                 fprintf(2,['Fail to remove. Folder ' path2file ' does not exist or has been removed\n']);
%             end
%         end
%         
%     end
% end

% Cleaning files related with design (unnecesary??)
tmpfiles = dir([filepath filesep 'LIMO_' filename]);
tmpfiles = {tmpfiles.name};

filecell = strfind(tmpfiles,['GLM' designname ]);
indxtmp = find(not(cellfun('isempty',filecell)));
if ~isempty(indxtmp)
    for nfile = 1:length(indxtmp)
        try
            file2delete = [filepath filesep 'LIMO_' filename filesep tmpfiles{indxtmp(nfile)}];
            if isdir(file2delete)
                rmdir(file2delete,'s');
            else
                delete(file2delete);
            end
        catch
            fprintf(2,['Fail to remove. File ' tmpfiles{indxtmp(nfile)} ' does not exist or has been removed\n']);
        end
    end
end
