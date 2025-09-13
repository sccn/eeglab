% POPHELP - Same as matlab HTHELP but does not crash under windows.
%
% Usage: >> pophelp( function );
%        >> pophelp( function, nonmatlab );
%
% Inputs:
%   function  - string for a Matlab function name 
%               (with or without the '.m' extension).
%   nonmatlab - [0|1], 1 the file is not a Matlab file
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: EEGLAB 

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function pophelp( funct, nonmatlab )

    if nargin <1
	    help pophelp;
	    return;
    end
    if nargin <2
	    nonmatlab = 0;
    end
    
    tmpvers = version;
    if str2double(tmpvers(1:2)) > 24 || ~isempty(findstr(tmpvers, 'R2024b'))
        showFunctionHelp(funct);
        return
    end
    
    if exist('help2html') == 2
        if length(funct) > 3 && strcmpi(funct(end-3:end), '.txt')
            %web(funct);
            fid = fopen(funct, 'r');
            text1 = textscan(fid, '%s', 'delimiter', '');
            fclose(fid);
            text1 = cellfun(@(x)[10 x], text1{1}, 'uniformoutput', false);
            tmp = char('text://<pre>', text1{:});
            tmp = tmp';
            tmp = tmp(:);
            web( tmp' );
        else
            pathHelpHTML = fileparts(which('help2html'));
            if ~isempty(findstr('NFT', pathHelpHTML)), rmpath(pathHelpHTML); end
            text1 = help2html(funct);
            if length(funct) > 4 && strcmpi(funct(1:4), 'pop_')
                try,
                    text2 = help2html(funct(5:end));
                    text1 = [text1 '<br><pre>___________________________________________________________________' 10 ...
                                   ' ' 10 ...
                                   ' The ''pop'' function above calls the eponymous Matlab function below' 10 ...
                                   ' and could use some of its optional parameters' 10 ...
                                   '___________________________________________________________________</pre><br><br>' text2 ];
                catch, end
            end
            if  exist('OCTAVE_VERSION','builtin') == 0
                web([ 'text://' text1 ]);
            else
                disp(text1);
            end
        end
    else
        [~,~,fileext] = fileparts(funct);
        if exist('OCTAVE_VERSION','builtin') == 0 && isequal(fileext, 'm')
            if usejava('jvm') && usejava('awt') % display available
                doc(funct);
                return;
	    end
        end
    
        if isempty(funct), return; end
        doc1 = readfunc(funct, nonmatlab);
        if length(funct) > 4 && strcmpi(funct(1:4), 'pop_')
            try
                doc2 = readfunc(funct(5:end), nonmatlab);
                doc1 = { doc1{:} ' ' ...
                    ' _________________________________________________________________ ' ...
                               ' ' ...
                               ' The ''pop'' function above calls the eponymous Matlab function below, ' ...
                               ' which may contain more information for some parameters. '...
                               ' ' ...
                               ' _________________________________________________________________ ' ...
                               ' ' ...
                        doc2{:} };
            catch, end
        end
    
        if exist('OCTAVE_VERSION','builtin') ~= 0
            for iRow = 1:length(doc1)
                disp(doc1{iRow});
            end
            return
        end
    
        % write file for help only
        icadefs;
        if ~isempty(EEGOPTION_PATH) % in icadefs above
            homefolder = EEGOPTION_PATH;
        elseif ispc
            homefolder = getenv('USERPROFILE');
        else homefolder = '~';
        end
            
        [~,funct] = fileparts(funct);
        fileTmp = fullfile(homefolder, [ funct '_doc.m' ]);
        fid = fopen(fileTmp, 'w');
        if fid ~= -1 
            for iDoc = 1:length(doc1)
                fprintf(fid, '%%%s\n', doc1{iDoc});
            end
            fclose(fid);
            doc(fileTmp)
            drawnow; pause(10);
            delete(fileTmp);
        else
            textgui(doc1);
            h = findobj('parent', gcf, 'style', 'slider');
            try 
                icadefs; 
            catch
                GUIBUTTONCOLOR = [0.8 0.8 0.8]; 
                GUITEXTCOLOR   = 'k'; 
            end
            set(h, 'backgroundcolor', GUIBUTTONCOLOR);
            h = findobj('parent', gcf, 'style', 'pushbutton');
            set(h, 'backgroundcolor', GUIBUTTONCOLOR);
            h = findobj('parent', gca);
            set(h, 'color', GUITEXTCOLOR);
            set(gcf, 'color', BACKCOLOR);
        end
    
    end
    return;
end

function [doc] = readfunc(funct, nonmatlab)
    
    doc = {};
    if iseeglabdeployed
        warndlg2([ 'Some help menus not available in compiled version.' 10 'Look up help online.' ] );
    end
    if nonmatlab	
	    fid = fopen( funct, 'r');
    else
	    if findstr( funct, '.m')
		    fid = fopen( funct, 'r');
	    else
		    fid = fopen( [funct '.m'], 'r');
	    end
    end
    
    if fid == -1
	    error('File not found');
    end
    
    sub = 1;
    try, 
        if ~isunix, sub = 0; end
    catch, end
    
    if nonmatlab
	    str = fgets( fid );
	    while ~feof(fid)
		    str = deblank(str(1:end-sub));
            doc = { doc{:} str(1:end) };    
            str = fgets( fid );
	    end
    else
	    str = fgets( fid );
	    while (str(1) == '%')
		    str = deblank(str(1:end-sub));
            doc = { doc{:} str(2:end) };    
		    str = fgets( fid );
	    end
    end
    fclose(fid);
end

function showFunctionHelp(fun, titleStr)
%SHOWFUNCTIONHELP  Render a function's HELP text in a modern popup.
%   showFunctionHelp(@fft)             % or showFunctionHelp('fft')
%   showFunctionHelp(@fft,'FFT Help')

    if nargin < 2 || isempty(titleStr)
        titleStr = 'Function Help';
    end
    if isa(fun,'function_handle')
        funName = func2str(fun);
    else
        funName = char(fun);
    end

    % Get HELP text exactly as MATLAB would show it
    txt = help(funName);
    if isempty(txt)
        txt = sprintf('No help found for "%s".', funName);
    end

    % Escape HTML special chars but keep newlines and spacing
    htmlBody = ['<pre>' localEscapeHTML(txt) '</pre>'];

    % Optional: prepend H1 line as a heading when present
    h1.file = which(funName);
    h1.text = upper(funName);
    if ~isempty(h1)
        htmlBody = sprintf('<h1>%s</h1><div class="h1file">%s</div>%s', ...
                           localEscapeHTML(strtrim(h1.text)), ...
                           localEscapeHTML(h1.file), htmlBody);
        titleStr = sprintf('%s — %s', funName, strtrim(h1.text)); %#ok<*NASGU>
    else
        htmlBody = sprintf('<h1>%s</h1>%s', localEscapeHTML(funName), htmlBody);
    end

    % Build minimal HTML shell with system fonts and dark mode
    html = [ ...
        "<!doctype html><html><head><meta charset=""utf-8"">" + ...
        "<meta name=""viewport"" content=""width=device-width, initial-scale=1"">" + ...
        "<style>" + ...
        "html,body{height:100%;} body{margin:16px;font-family:system-ui,-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;line-height:1.4;font-size:13px;}" + ...
        "h1{margin:0 0 8px 0;font-size:16px;}" + ...
        ".h1file{color:#666;margin:0 0 12px 0;font-size:12px;}" + ...
        "pre{white-space:pre-wrap;word-wrap:break-word;margin:0;font-family:Consolas,Menlo,Monaco,monospace;font-size:12px;}" + ...
        "@media (prefers-color-scheme: dark) { body{background:#121212;color:#E0E0E0;} .h1file{color:#9AA0A6;} }" + ...
        "</style></head><body>" + htmlBody + "</body></html>" ...
        ];

    % Create popup window
    f = uifigure('Name', titleStr, 'Position', [100 100 720 520], ...
                 'AutoResizeChildren','on');
    g = uigridlayout(f, [1 1]); g.RowHeight = {"1x"}; g.ColumnWidth = {"1x"};
    h = uihtml(g);
    h.HTMLSource = char(html); % uihtml accepts char or string

    % Local helpers
    function s = localEscapeHTML(s)
        s = strrep(s,'&','&amp;');
        s = strrep(s,'<','&lt;');
        s = strrep(s,'>','&gt;');
    end

    function h1 = getH1Line(name)
        h1 = struct('text','','file','');
        try
            fpath = which(name);
            if isempty(fpath), return, end
            fid = fopen(fpath,'r');
            if fid<0, return, end
            cleaner = onCleanup(@() fclose(fid));
            line = fgetl(fid);
            % MATLAB-style H1 is on the first comment line after the function name
            % pattern: % H1 summary
            while ischar(line)
                t = strtrim(line);
                if startsWith(t,'%')
                    t = strtrim(erase(t,'%'));
                    if ~isempty(t)
                        h1.text = t;
                        h1.file = fpath;
                        return
                    end
                elseif ~isempty(t)
                    return
                end
                line = fgetl(fid);
            end
        catch
            % best effort only
        end
    end
end