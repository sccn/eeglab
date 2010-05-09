% plotchans3d() -  Plots the 3-D configuration from a Polhemus ELP
%                  file. The axes of the Cartesian coordinate system are
%                  plotted. The nose is plotted as a bold red segment.
% Usage:
%        >> plotchans3d( elpfile, zs);
%        >> plotchans3d( [X,Y,Z], elecnames, zs);
%
% Inputs:
%        elpfile - Polhemus ELP (electrode position) file 
%        [X,Y,Z] - array of 3-D coordinates
%        elecnames - cell array containing the names of the electrodes
%
% Optional input:
%        zs - vector of electrode indices giving the subset to print labels for 
%             in the plot
%
% Author: Luca A. Finelli, SCCN/INC/UCSD, La Jolla, 02/2002
%  
% See also: topoplot(), readlocs(), readelp(), 
%           polhemus2topo(), pop_chanedit()
        
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

% 04/01/02 changed header for help2html compatibility -ad
% 04/01/02 debuging zs -ad
% 04/01/02 adding extra optional [X,Y,Z] argument -ad

function  plotchans3d(elpfile, arg2, arg3)

if nargin<1
    help plotchans3d;
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
    if nargin > 1
        elocname = arg2;
    else 
        elocname = [];
    end;
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

if ~isempty(elocname)
    plot3(X(zs),Y(zs),Z(zs),'b*')
end;

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
if ~isempty(elocname)
    for i = 1:length(zs)
        text(X(zs(i)),Y(zs(i)),Z(zs(i))+eps,elocname(zs(i)),'HorizontalAlignment','center',...
             'VerticalAlignment','middle','Color',[0 0 0],...
             'FontSize',10)
    end
end;
%axis(repmat([-lim lim],1,3))
axis([-lim lim -lim lim -lim*0.5 lim])
axis equal;
rotate3d on
try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end
