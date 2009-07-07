% XMLSTRUCT - Reads XML data from file and list as structures.
%
% st=xmlstruct('xmlstring');
%
%  This returns a long array of struct describing all XML tags and
%  and text between the tags. All parameters belonging to the tag
%  are stored in an a sub structure '___.data'. 
%
%
% st=xmlstruct('xmlstring','sub');
%
%  This organizes all structure in a tree shape with sub structures.
%  Then a opening tag (type=1) opens up a sub level and puts them
%  all following data in '__.sub' and the corresponding closing
%  tag exits that sub level.
%  
%
% SMALL EXAMPLE 1:
%
%   st=xmlstruct('<creator type="composer">Franz Schubert</creator>')
%
%-st
% |-st(1)
% | |-name .............. creator
% | |-data
% | | \-type .............. composer
% | \-type ..............      1
% |-st(2)
% | |-name .............. []
% | |-data .............. Franz Schubert
% | \-type ..............      0
% \-st(3)
%   |-name .............. creator
%   |-data .............. []
%   \-type ..............      3
%
% SMALL EXAMPLE 2:
%   
%   st=xmlstruct('<creator type="composer">Franz Schubert</creator>','sub')
%
%-st
% |-name .............. creator
% |-data
% | \-type .............. composer
% |-type ..............      1
% \-sub
%   |-sub(1)
%   | |-name .............. []
%   | |-data .............. Franz Schubert
%   | \-type ..............      0
%   \-sub(2)
%     |-name .............. creator
%     |-data .............. []
%     \-type ..............      3
%     
% Tip, You can quick and easy read a XML-file into a matlab string with:
%
%   f=fopen('myfile.xml');
%   xmlstring=char(fread(f,[1 inf],'uint8'));
%   fclose(f);
%   
% SE ALSO: __XMLDATA__
%
% Look at:  http://www.rydesater.com/tools OR http://petrydpc.ite.mh.se/tools
%
%
%   LICENSE:
%  ========= 
%
%   Copyright (C) Peter Rydesäter 2002, Mitthögskolan, SWEDEN
%
%   This program is free software; you can redistribute it and/or
%   modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation; either version 2
%   of the License, or (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
%   
%
%   Please, send me an e-mail if you use this and/or improve it. 
%   If you think this is usful fore your work, please remember me and my
%   hard work with it!
%
%   Peter.Rydesater@mh.se
%
%   Peter Rydesäter
%   Mitthögskolan (Mid Sweden University)
%   SE-831 25 ÖSTERSUND
%   SWEDEN
%


function [st,p]=xmlstruct(str,par,p)
    if nargin<3,
      p=1;
    end
    n=0;
    while p,
      n=n+1;
      if mod(n,500)==1,
	st(n+501).name='';
      end      
      [a,b,c,p]=xmldata(str,p);    
      st(n).name=a;
      st(n).data=b;
      st(n).type=c;      
      if nargin>1,
	if c==1,
	  [st(n).sub,p]=xmlstruct(str,'sub',p);
	elseif c==3,
	  if nargin==3, break; end	  
	end
      end
    end
    if ischar(st(n).name),
      st=st(1:n);
    else
      st=st(1:n-1);
    end
    return;
    
