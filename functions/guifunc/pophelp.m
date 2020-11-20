% pophelp() - Same as matlab HTHELP but does not crash under windows.
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
% See also: eeglab() 

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

function pophelp( funct, nonmatlab );

if nargin <1
	help pophelp;
	return;
end
if nargin <2
	nonmatlab = 0;
end

if exist('help2html')
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

        web([ 'text://' text1 ]);
    end
else
    if isempty(funct), return; end
    doc1 = readfunc(funct, nonmatlab);
    if length(funct) > 4 && strcmpi(funct(1:4), 'pop_')
        try,
            doc2 = readfunc(funct(5:end), nonmatlab);
            doc1 = { doc1{:} ' _________________________________________________________________ ' ...
                           ' ' ...
                           ' The ''pop'' function above calls the eponymous Matlab function below, ' ...
                           ' which may contain more information for some parameters. '...
                           ' ' ...
                           ' _________________________________________________________________ ' ...
                           ' ' ...
                    doc2{:} };
        catch, end
    end

    textgui(doc1);1000
    h = findobj('parent', gcf, 'style', 'slider');
    try, icadefs; catch, 
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
return;

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
