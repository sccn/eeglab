function scptest(p,fid)
% SCPTEST checks files for SCP-ECG dataformat. 
%
%  scptest(p)
%	checks recursively all files in the subdirectory p
%  scptest(f)
%	checks file f 
%
%  fid = fopen('logfile.log','w');
%  scptest(p,fid)
%     redirects output into logfile with open handle fid. 
%
%
% see also: SOPEN, SCLOSE, SCPOPEN
%
% References: 
% [1] The OpenECG project: http://www.openecg.org/ 
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

%	$Id: scptest.m,v 1.1 2009-01-30 06:04:39 arno Exp $
%	Copyright (C) 2004 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

if nargin<2, fid=1; end;

fn = dir(fullfile(p,'*')); 
for k = 3:length(fn),
        f = fullfile(p,fn(k).name);
        if fn(k).isdir
                scptest(f,fid);
        else
		FLAG = 0; 
                try,
                        tic,
                        H = sopen(f);
                        H = sclose(H);
			%drawnow
			FLAG = strcmp(H.TYPE,'SCP');
		catch;

		end;

		if FLAG,
                        fprintf(fid,'%5.2f\tYES\t%s\n',toc,f);
                else
                        fprintf(fid,'%5.2f\t---\t%s\n',toc,f);
                end;
        end;
end;
