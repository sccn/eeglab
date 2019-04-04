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
