% testelp() -  Plots the 3D configuration from a Polhemus ELP
%           file. The axis of the cartesian coordinate system are
%           plotted. The nose is plotted as a bold red segment.
%
% Usage:
%        >> testelp( elpfile, zs)
%        >> testelp( [X,Y,Z], elecnames, zs)
%
% Inputs:
%        elpfile = ELP file from Polhemus.
%        [X,Y,Z] = array of 3D coordinates
%        elecnames = cell array containing the name of the electrodes
%
% Optional input:
%        zs = vector of electrode subset to print labels
%
% Author: Luca A. Finelli, SCCN/INC/UCSD, La Jolla, 02/2002
%  
% See also:
%        topoplot(), readlocs(), readelp(), polhemus2topo()
        
%123456789012345678901234567890123456789012345678901234567890123456789012

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: not supported by cvs2svn $
% Revision 1.1  2002/04/05 17:48:29  jorn
% Initial revision
%

% 04/01/02 changed header for help2html compatibility -ad
% 04/01/02 debuging zs -ad
% 04/01/02 adding extra optional [X,Y,Z] argument -ad

function  testelp(elpfile, arg2, arg3)

if nargin<1
    help testelp;
    return;
end;
zs = [];
if isstr(elpfile)
    if nargin > 1
        zs = arg2;
    end;
    [elocstruct, elocname, X, Y, Z ] =readelp([elpfile]);
else
    X = elpfile(:,1)';
    Y = elpfile(:,2)';
    Z = elpfile(:,3)';
    elocname = arg2;
    if nargin > 2
        zs = arg3;
    end;
end;

if isempty(zs)
   zs = [1:length(elocname)];
end;

%zs =[3 7 15 26 36 46 56 64 69 71 72];

figure
lim=1.05*max([X Y Z]);
eps=lim/20;
plot3(X,Y,Z,'ro')
hold on

elocname(zs);
plot3(X(zs),Y(zs),Z(zs),'b*')

plot3([0.08 0.12],[0 0],[0 0],'r','LineWidth',4) % nose
plot3([0 lim],[0 0],[0 0],'b--')                 % axis
plot3([0 0],[0 lim],[0 0],'g--')
plot3([0 0],[0 0],[0 lim],'r--')
plot3(0,0,0,'b+')
text(lim+eps,0,0,'X','HorizontalAlignment','center',...
	'VerticalAlignment','middle','Color',[0 0 0],...
	'FontSize',10)
text(0,lim+eps,0,'Y','HorizontalAlignment','center',...
	'VerticalAlignment','middle','Color',[0 0 0],...
	'FontSize',10)
text(0,0,lim+eps,'Z','HorizontalAlignment','center',...
	'VerticalAlignment','middle','Color',[0 0 0],...
	'FontSize',10)
box on
for i = 1:length(zs)
    text(X(zs(i)),Y(zs(i)),Z(zs(i))+eps,elocname(zs(i)),'HorizontalAlignment','center',...
	'VerticalAlignment','middle','Color',[0 0 0],...
	'FontSize',10)
end
%axis(repmat([-lim lim],1,3))
axis([-lim lim -lim lim -lim*0.5 lim])
rotate3d on
