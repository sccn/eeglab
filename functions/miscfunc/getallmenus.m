% getallmenus() - get all submenus of a window or a menu and return 
%                 a tree.
%
% Usage:
%   >> [tree nb] = getallmenus( handler );
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

function [txt, nb, labels] = getallmenus( handler, level )

NBBLANK = 6; % number of blank for each submenu input

	if nargin < 1
		help getallmenus;
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
		for i=1:length(	allmenu );
			[txtmp nbtmp tmplab] = getallmenus(allmenu(i), level+1);
			txtmp = [ blanks(level*NBBLANK) txtmp ];
			txt = [ txtmp txt ];
			labels = { tmplab labels{:} };
			nb = nb+nbtmp;
		end
	end
	try
		txt = [ get(handler, 'Label') 10 txt ];
		nb = nb+1;
	catch, end
	if isempty(labels)
		labels = { nb };
	end;	

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
return;
