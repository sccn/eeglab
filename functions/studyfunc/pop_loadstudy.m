% pop_loadstudy() - load an existing STUDY set and its corresponding ALLEEG structure
%
% Usage:
%   >> [STUDY ALLEEG] = pop_loadstudy; % use an interactive pop-up window 
%   >> [STUDY ALLEEG] = pop_savestudy( 'key', 'val', ...); % no pop-up
%                                              
%
% Optional inputs:
%   'filename' - [string] name of the study set file to load
%   'filepath' - [string] path of the study set file to load
%
% Outputs:
% STUDY - an existing STUDY set loaded.
% ALLEEG - an ALLEEG structure with the datasets corresponding to the
%                datasets in the study set.    
%
function [STUDY, ALLEEG] = pop_loadstudy(varargin)

if isempty(varargin)
    [filename, filepath] = uigetfile2('*.study', 'Load a STUDY -- pop_loadstudy()'); 
    if ~strncmp(filename(end-5:end), '.study',6)
        if isempty(strfind(filename,'.'))
            filename = [filename '.study'];
        else
            filename = [filename(1:strfind(filename,'.')-1) '.study'];
        end
    end
else
    for k = 1:2:length(varargin)
        switch varargin{k}
            case 'filename'
                filename = varargin{k+1};
            case 'filepath'
                filepath = varargin{k+1};
        end
    end
end

if (~isempty(filename)) & (~isempty(filepath))
    STUDYfile = fullfile(filepath,filename);
    try 
        eval(['load ' STUDYfile ' -mat']);
    catch
        error(['pop_loadstudy: the Study file ' STUDYfile ' could not be loaded, check file name and file path']);
    end
else
    error(['pop_loadstudy: No Study file to load was provided.']);
end
  
ALLEEG = load_ALLEEG(STUDY);

    
       