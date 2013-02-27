% makehelpfiles() - generate function help pages 
%
% Usage: 
%   >> makehelpfiles(list);
%
% Input:
%   'folder'     - [string] folder name to process
%   'outputfile' - [string] file in which to write the help
%   'title'      - [string] title for the help
%
% Author: Arnaud Delorme, UCSD, 2013
%
% Example: Generate EEGLAB help menus for adminfunc folder
%  makehelpfiles('folder', 'adminfunc'   ,'title', 'Admin functions', 'outputfile','adminfunc/eeg_helpadmin.m' );          
%  makehelpfiles('folder', 'guifunc'     ,'title', 'Graphic interface builder functions', 'outputfile','adminfunc/eeg_helpgui.m' );          
%  makehelpfiles('folder', 'miscfunc'    ,'title', 'Miscelaneous functions not used by EEGLAB graphic interface', 'outputfile','adminfunc/eeg_helpmisc.m' );          
%  makehelpfiles('folder', 'popfunc'     ,'title', 'EEGLAB graphic interface functions', 'outputfile','adminfunc/eeg_helppop.m' );          
%  makehelpfiles('folder', 'sigprocfunc' ,'title', 'EEGLAB signal processing functions', 'outputfile','adminfunc/eeg_helpsigproc.m' );         
%  makehelpfiles('folder', 'statistics'  ,'title', 'EEGLAB statistics functions', 'outputfile','adminfunc/eeg_helpstatistics.m' );         
%  makehelpfiles('folder', 'studyfunc'   ,'title', 'EEGLAB group processing (STUDY) functions', 'outputfile','adminfunc/eeg_helpstudy.m' );         
%  makehelpfiles('folder', 'timefreqfunc','title', 'EEGLAB time-frequency functions', 'outputfile','adminfunc/eeg_helptimefreq.m' );         

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 2002
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

function makehelpfiles( varargin );

if nargin < 1
    help makehelpfiles;
    return;
end;
    
opt = finputcheck( varargin,  { 'folder'     'string'   { }  '';
                                'outputfile' 'string'   { }  '';
                                'title'      'string'   { }  '' }, 'makehelpfiles');
if isstr(opt), error(opt); end;
if isempty(opt.folder),     error('You need to specify a folder'); end;
if isempty(opt.outputfile), error('You need to specify an output file'); end;

fo = fopen( opt.outputfile, 'w');
if ~isempty(opt.title)
     fprintf(fo, '%%%s (%s folder):\n', opt.title, opt.folder);
else fprintf(fo, '%% *Content of %s folder:*\n', opt.folder);
end;
dirContent = dir(fullfile(opt.folder, '*.m'));
dirContent = { dirContent.name };
for iFile = 1:length(dirContent)
    fileName  = dirContent{iFile};
    fidTmp    = fopen(fullfile(opt.folder, fileName), 'r');
    firstLine = fgetl(fidTmp);
    
    % get help from the first line
    if isempty(firstLine) || firstLine(1) ~= '%'
        firstLine = fgetl(fidTmp);
    end;
    fclose(fidTmp);
        
    if isempty(firstLine) || firstLine(1) ~= '%'
        firstLineText = 'No help information'; 
    else
        indexMinus    = find(firstLine == '-');
        if ~isempty(indexMinus)
             firstLineText = deblank(firstLine(indexMinus(1)+1:end));
        else firstLineText = deblank(firstLine(2:end));
        end;
        if isempty(firstLineText), firstLineText = 'No help information'; end;
        if firstLineText(1) == ' ', firstLineText(1) = []; end;
        if firstLineText(1) == ' ', firstLineText(1) = []; end;
        if firstLineText(1) == ' ', firstLineText(1) = []; end;
        firstLineText(1) = upper(firstLineText(1));
        if firstLineText(end) ~= '.', firstLineText = [ firstLineText '...' ]; end;
    end;
    
    refFunction = sprintf('<a href="matlab:helpwin %s">%s</a>', fileName(1:end-2), fileName(1:end-2));
    fprintf(fo, '%%  %-*s - %s\n', 50+length(fileName(1:end-2)), refFunction, firstLineText);
end;
fclose( fo );
