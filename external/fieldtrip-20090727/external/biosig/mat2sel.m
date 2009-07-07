function HDR = mat2sel(FileName,fnout);
% MAT2SEL is useful for converting visual artifact selection 
%   with gBSanalyze into BKR/EPS compliant *.SEL files. 
%
%   HDR = mat2sel(infile, [selfile]);     
%      infile   input filename (*.BKR or *.MAT are supported)
%      selfile  output file with info artifact selection [OPTIONAL]
%               if no selfile is provided, the information will be 
%               stored with extension '.SEL' and the same filename and
%               path. 
%
%   The artifact selection is available in HDR.ArtifactSelection and 
%   can be also loaded with 
%       [s,HDR]=sload(filename);
%   Currently, artifact selection of  EPS/BKR-Software and 
%   gBSanalyze is supported. 
%  
% see also: SLOAD

%	$Revision: 1.1 $
%	$Id: mat2sel.m,v 1.1 2009-07-07 02:23:46 arno Exp $
%	Copyright (C) 2004 by Alois Schloegl 
%	a.schloegl@ieee.org	

% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 2 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.

[s,HDR] = sload(FileName);

if isfield(HDR,'ArtifactSelection'),
        artifact = HDR.ArtifactSelection;
        
        save([HDR.FILE.Name, '_artifact'], 'artifact');
        
        if nargin<2,
                fnout = [HDR.FILE.Name, '.sel']
                fnout = fullfile(HDR.FILE.Path,[HDR.FILE.Name, '.sel']);
        end
        if 0, exist(fnout)==2,
                fprintf(1,'File %s exists! Do you want to overwrite ',fnout);
                answer = input('[Y/N] ?');
        else
                answer = 'Y';
        end;
        if any(answer=='jJyYzZ'),
                fid = fopen(fnout,'w');
                fprintf(fid, '%i\r\n', artifact);
                fclose(fid);
        end;
end;