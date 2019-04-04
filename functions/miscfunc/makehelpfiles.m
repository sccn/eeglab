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

function makehelpfiles( varargin );

if nargin < 1
    help makehelpfiles;
    return;
end
    
opt = finputcheck( varargin,  { 'folder'     'string'   { }  '';
                                'outputfile' 'string'   { }  '';
                                'title'      'string'   { }  '' }, 'makehelpfiles');
if ischar(opt), error(opt); end
if isempty(opt.folder),     error('You need to specify a folder'); end
if isempty(opt.outputfile), error('You need to specify an output file'); end

fo = fopen( opt.outputfile, 'w');
if ~isempty(opt.title)
     fprintf(fo, '%%%s (%s folder):\n', opt.title, opt.folder);
else fprintf(fo, '%% *Content of %s folder:*\n', opt.folder);
end
dirContent = dir(fullfile(opt.folder, '*.m'));
dirContent = { dirContent.name };
for iFile = 1:length(dirContent)
    fileName  = dirContent{iFile};
    fidTmp    = fopen(fullfile(opt.folder, fileName), 'r');
    firstLine = fgetl(fidTmp);
    
    % get help from the first line
    if isempty(firstLine) || firstLine(1) ~= '%'
        firstLine = fgetl(fidTmp);
    end
    fclose(fidTmp);
        
    if isempty(firstLine) || firstLine(1) ~= '%'
        firstLineText = 'No help information'; 
    else
        indexMinus    = find(firstLine == '-');
        if ~isempty(indexMinus)
             firstLineText = deblank(firstLine(indexMinus(1)+1:end));
        else firstLineText = deblank(firstLine(2:end));
        end
        if isempty(firstLineText), firstLineText = 'No help information'; end
        if firstLineText(1) == ' ', firstLineText(1) = []; end
        if firstLineText(1) == ' ', firstLineText(1) = []; end
        if firstLineText(1) == ' ', firstLineText(1) = []; end
        firstLineText(1) = upper(firstLineText(1));
        if firstLineText(end) ~= '.', firstLineText = [ firstLineText '...' ]; end
    end
    
    refFunction = sprintf('<a href="matlab:helpwin %s">%s</a>', fileName(1:end-2), fileName(1:end-2));
    fprintf(fo, '%%  %-*s - %s\n', 50+length(fileName(1:end-2)), refFunction, firstLineText);
end
fclose( fo );
