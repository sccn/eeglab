function [EVENT,cc] = adb2event(fn,Fs)
% ADB2EVENT loads and artifact scoring file of the 
%   artifact database (ADB) of sleep EEG [1]
%   and converts the data into event information 
% 
%  [s,H] = sload(filename.edf)     
%  H.EVENT = adb2event([H.FILE.Name,'.txt], H.SampleRate);
%
% see also: SOPEN, SREAD, SSEEK, STELL, SCLOSE, SWRITE, SEOF
%
% Reference(s):
% [1] Artifact database of sleep EEG. Available online http://www.dpmi.tu-graz.ac.at/ADB/
% [2] A. Schlögl, P. Anderer, M.-J. Barbanoj, G. Klösch,G. Gruber, J.L. Lorenzo, O. Filz, M. Koivuluoma, I. Rezek, S.J. Roberts,A. Värri, P. Rappelsberger, G. Pfurtscheller, G. Dorffner
%       Artifact processing of the sleep EEG in the "SIESTA"-project,
%       Proceedings EMBEC'99, Part II, pp.1644-1645, 4-7. Nov. 1999,Vienna, Austria.
% [3] A. Schlögl, P. Anderer, S.J. Roberts, M. Pregenzer, G.Pfurtscheller.
%       Artefact detection in sleep EEG by the use of Kalman filtering.
%       Proceedings EMBEC'99, Part II, pp.1648-1649, 4-7. Nov. 1999,Vienna, Austria.
% [4] A. Schlögl, P. Anderer, M.-J. Barbanoj, G. Dorffner, G. Gruber, G. Klösch, J.L. Lorenzo, P. Rappelsberger, G. Pfurtscheller.
%       Artifacts in the sleep EEG - A database for the evaluation of automated processing methods.
%       Proceedings of the Third International Congress of the World Federation of Sleep Research Societies (WFSRS). Editors: H. Schulz. P.L. Parmeggiani, and M. Chase. Sleep Research Online 1999:2 (Supplement 1), p. 586. 
%       available online: http://www.sro.org/cftemplate/wfsrscongress/indiv.cfm?ID=19998586


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
%	$Id: adb2event.m,v 1.1 2009-07-07 02:23:45 arno Exp $
%	(C) 1997-2004 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

ma = load(fn);

% code from MAK2BIN.M (C) 1998-2004 A. Schlögl 
ERG = zeros(size(ma));

%%%% one artifact %%%%
for k=0:9,
        if exist('OCTAVE_VERSION')==5
                ERG = ERG+(ma==k)*2^k;
        else
                ERG(ma==k) = 2^k;
        end;
end;

%%%% more than one artifact %%%%
[i,j] = find(ma>9);
L='123456789';
for k=1:length(i),
        b=int2str(ma(i(k),j(k)));
        erg=0;
        for l=1:9,
                if any(b==L(l)), erg=erg+2^l; end;
        end;        
        ERG(i(k),j(k)) = erg;
end;

N   = 0;
POS = [];
TYP = [];
DUR = [];
CHN = [];
cc  = zeros(1,10);
for k = 1:9,
        for c = 1:7;%size(ERG,2),
                tmp = [0;~~(bitand(ERG(:,c),2^k));0];
 
                cc(k+1) = cc(k+1) + sum(tmp);
                pos = find(diff(tmp)>0);
                pos2 = find(diff(tmp)<0);
                n   = length(pos);
                
                POS = [POS; pos(:)];
                TYP = [TYP; repmat(k,n,1)];
                CHN = [CHN; repmat(c,n,1)];
                DUR = [DUR; pos2(:)-pos(:)];
                N   = N + n;
        end;
end;

EVENT.Fs = 1;
if nargin>1,
        EVENT.Fs = Fs;
end;

[tmp,ix] = sort(POS);
EVENT.POS = (POS(ix)-1)*EVENT.Fs+1;
EVENT.TYP = TYP(ix) + hex2dec('0100');
EVENT.CHN = CHN(ix);
EVENT.DUR = DUR(ix)*EVENT.Fs;
EVENT.N   = N;

%EVENT.ERG = ERG; 