% pop_savestudy() - save a STUDY structures
%
% Usage:
%   >> STUDY = pop_savestudy( STUDY ); % use an interactive pop-up window 
%   >> STUDY = pop_savestudy( STUDY, 'key', 'val', ...); % no pop-up
%                                              
% Inputs:
%   STUDY       - STUDY structure. 
%
% Optional inputs:
%   'filename' - [string] name of the file to save to {default: STUDY.filename}
%   'filepath' - [string] path of the file to save to {default: STUDY.filepath}
% 
function [STUDY] = pop_savestudy(STUDY, varargin);

if nargin < 1
	help pop_savestudy;
	return;
end;
if isempty(STUDY)  , error('Cannot save empty STUDY'); end;
if length(STUDY) >1, error('This function does not save multiple STUDY sets'); end;

if nargin < 2
    % pop up window to ask for file type
    % ----------------------------------
    [g.filename, g.filepath] = uiputfile2('*.study', 'Save STUDY with .study extension -- pop_savestudy()'); 
    if ~strncmp(g.filename(end-5:end), '.study',6)
        if isempty(strfind(g.filename,'.'))
            g.filename = [g.filename '.study'];
        else
            g.filename = [g.filename(1:strfind(g.filename,'.')-1) '.study'];
        end
    end
end

if nargin > 2 
    options = varargin;
    g = finputcheck(options,  { 'filename'   'string'   []     STUDY.filename;
                            'filepath'   'string'   []     STUDY.filepath;});
    if isstr(g), error(g); end;
end
if isempty(g.filename)
    disp('pop_savestudy: no filename is available please make sure there is a filename');
    return;
end
if ~strncmp(g.filename(end-5:end), '.study',6)
    if isempty(strfind(g.filename,'.'))
        g.filename = [g.filename '.study'];
    else
        g.filename = [g.filename(1:strfind(g.filename,'.')-1) '.study'];
    end
end

%    [filepath filenamenoext ext] = fileparts(varargin{1});
%    filename = [filenamenoext '.study']; % make sure a .study extension
    
STUDY.filepath = g.filepath;
STUDY.filename = g.filename;
STUDYfile = fullfile(STUDY.filepath,STUDY.filename);
try 
    STUDYfile = fullfile(STUDY.filepath,STUDY.filename);
    eval(['save ' STUDYfile ' STUDY -mat']);
catch
    try
        eval(['save ' STUDY.filepath '/' STUDY.filename ' STUDY -mat']);
        STUDY.filepath = [ STUDY.filepath '/' ];
    catch
        try 
            eval(['save ' STUDY.filepath '\' STUDY.filename ' STUDY -mat']);
            STUDY.filepath = [ STUDY.filepath '\' ];
        catch
            warning('pop_preclust: STUDY was not saved check file path and file name')
        end
    end
end