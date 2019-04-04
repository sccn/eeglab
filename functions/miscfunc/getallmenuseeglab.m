% getallmenuseeglab() - get all submenus of a window or a menu and return 
%                 a tree. The function will also look for callback.
%
% Usage:
%   >> [tree nb] = getallmenuseeglab( handler );
%
% Inputs:
%   handler    - handler of the window or of a menu
%
% Outputs:
%   tree       - text output
%   nb         - number of elements in the tree
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001

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

function [txt, nb, labels] = getallmenuseeglab( handler, level )

    NBBLANK = 6; % number of blank for each submenu input

	if nargin < 1
		help getallmenuseeglab;
		return;
	end
	if nargin < 2
		level = 0;
	end
	    
	txt = '';
	nb = 0;
	labels = {};
	allmenu = findobj('parent', handler, 'type', 'uimenu');
	allmenu = allmenu(end:-1:1);
	if ~isempty(allmenu)
		for i=length(allmenu):-1:1
			[txtmp nbtmp tmplab] = getallmenuseeglab(allmenu(i), level+1);
			txtmp = [ '% ' blanks(level*NBBLANK) txtmp ];
			txt = [ txtmp txt ];
			labels = { tmplab labels{:} };
			nb = nb+nbtmp;
		end
	end
	try
        lab = get(handler, 'Label');
        cb  = get(handler, 'Callback');
        cb  = extract_callback(cb);
        if ~isempty(cb)
    		 newtxt = [ lab ' - <a href="matlab:helpwin ' cb '">' cb '</a>'];
        else newtxt = [ lab ];
        end
        txt = [ newtxt 10  txt ];
		%txt = [ get(handler, 'Label') 10 txt ];
		nb = nb+1;
	catch, end
	if isempty(labels)
		labels = { nb };
	end;	
    if level == 0
        fid = fopen('tmpfile.m', 'w');
        fprintf(fid, '%s', txt);
        fclose(fid);
        disp(' ');
        disp('Results saved in tmpfile.m');
    end

% transform into array of text
% ----------------------------
if nargin < 2
	txt = [10 txt];
	newtext = zeros(1000, 1000);
	maxlength = 0;
	lines = find( txt == 10 );
	for index = 1:length(lines)-1
		tmptext = txt(lines(index)+1:lines(index+1)-1); 	
		if maxlength < length( tmptext ), maxlength = length( tmptext ); end
		newtext(index, 1:length(tmptext)) = tmptext;
	end
	txt = char( newtext(1:index+1, 1:maxlength) );
end;		

% extract plugin name
% -------------------
function cbout = extract_callback(cbin);

funcList = { 'pop_' 'topoplot' 'eeg_' };
for iList = 1:3
    indList = findstr(funcList{iList}, cbin);
    if ~isempty(indList), 
        if strcmpi(cbin(indList(1):indList(1)+length('pop_stdwarn')-1), 'pop_stdwarn')
            indList = findstr(funcList{iList}, cbin(indList(1)+1:end))+indList(1);
        end
        break; 
    end
end
if ~isempty(indList)
    indEndList = find( cbin(indList(1):end) == '(' );
    if isempty(indEndList) || indEndList(1) > 25
        indEndList = find( cbin(indList(1):end) == ';' );
        if cbin(indList(1)+indEndList(1)-2) == ')'
            indEndList = indEndList-2;
        end
    end
    cbout = cbin(indList(1):indList(1)+indEndList(1)-2);
else
    cbout = '';
end
