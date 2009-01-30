function v = biosigVersion;
% BIOSIGVERSION returns the version number of 
%  the BioSig-toolbox
%

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

%	$Id: biosigVersion.m,v 1.1 2009-01-30 06:04:40 arno Exp $
%	Copyright 2007 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/


global BIOSIG_GLOBAL;
if ~isfield(BIOSIG_GLOBAL,'VERSION')
	BIOSIG_GLOBAL.VERSION = 0 ; 
end; 
if ~BIOSIG_GLOBAL.VERSION,
	p = which('biosigVersion');
	[p1,f1,e1] = fileparts(p);
	[p2,f2,e2] = fileparts(p1);
	fn = dir(fullfile(p2,'VERSION')); 
	fid = fopen(fullfile(p2,'VERSION'),'rt'); 
	while ~feof(fid),
		line = fgetl(fid); 
		if strncmp(line,'# Version',9); 
			[t,r]=strtok(line,[9,':']);
			[t,r]=strtok(r,[9,':']);
			BIOSIG_GLOBAL.VERSION = t; 
		end; 	
	end; 	
	fclose(fid); 
end; 
v = BIOSIG_GLOBAL.VERSION; 

