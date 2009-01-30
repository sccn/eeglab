function [EVENT,cc] = wscore2event(f0,fc)
% WSCORE2EVENT loads WSCORE event files 
%   and converts the data into event information 
%
% EVENT = wscore2event(f0,fc)
%       f0 event file 
%       fc event codes
%
% see also: SOPEN, SREAD, SSEEK, STELL, SCLOSE, SWRITE, SEOF
%
% Reference(s):
% [1] Artifact database of sleep EEG. Available online http://www.dpmi.tu-graz.ac.at/ADB/
%
% [2] A. Schlögl, P. Anderer, M.-J. Barbanoj, G. Klösch,G. Gruber, J.L. Lorenzo, O. Filz, M. Koivuluoma, I. Rezek, S.J. Roberts,A. Värri, P. Rappelsberger, G. Pfurtscheller, G. Dorffner
%       Artifact processing of the sleep EEG in the "SIESTA"-project,
%       Proceedings EMBEC'99, Part II, pp.1644-1645, 4-7. Nov. 1999,Vienna, Austria.
% [3] A. Schlögl, P. Anderer, S.J. Roberts, M. Pregenzer, G.Pfurtscheller.
%       Artefact detection in sleep EEG by the use of Kalman filtering.
%       Proceedings EMBEC'99, Part II, pp.1648-1649, 4-7. Nov. 1999,Vienna, Austria.
% [4] A. Schlögl, P. Anderer, M.-J. Barbanoj, G. Dorffner, G. Gruber, G. Klösch, J.L. Lorenzo, P. Rappelsberger, G. Pfurtscheller.
%       Artifacts in the sleep EEG - A database for the evaluation of automatedprocessing methods.
%       Proceedings of the Third International Congress of the World Federation of Sleep Research Societies (WFSRS). Editors: H. Schulz. P.L. Parmeggiani, and M. Chase. Sleep Research Online 1999:2 (Supplement 1), p. 586. 



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

%	$Revision: 1.1 $
%	$Id: wscore2event.m,v 1.1 2009-01-30 06:04:42 arno Exp $
%	(C) 1997-2004 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/


if nargin>1,
        k = 0;
        fid = fopen(fc,'r')
        while ~feof(fid),
                k = k+1;
                s=fgetl(fid);
                [t1,s]=strtok(s);
                %[t2,s]=strtok(s);
                EVENT.CodeDesc{k,1}=s(2:end);
                EVENT.CodeIndex(k,1)=str2double(t1);
        end;
        fclose(fid);
end;

k = 0;
fid = fopen(f0,'r');
while ~feof(fid),
      k = k+1;
      s=fgetl(fid);
      [t1,s]=strtok(s);
      [t2,s]=strtok(s);
      EVENT.POS(k,1)=str2double(t1);
      EVENT.TYP(k,1)=str2double(t2);
end;
fclose(fid);

