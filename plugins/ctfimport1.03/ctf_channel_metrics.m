function [stats] = ctf_channel_metrics(ctf)

% ctf_channel_metrics - CTF channel geometry
% 
% [stats] = ctf_channel_metrics(ctf)
%
% stats output contains measures of gradiometer separation and
% inter-gradiometer distances.
%
%      <>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%      <                                                      >
%      <                    DISCLAIMER:                       >
%      <                                                      >
%      < THIS PROGRAM IS INTENDED FOR RESEARCH PURPOSES ONLY. >
%      < THIS PROGRAM IS IN NO WAY INTENDED FOR CLINICAL OR   >
%      <                    OFFICIAL USE.                     >
%      <                                                      >
%      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>
%


% $Revision: 1.1 $ $Date: 2009-01-30 03:49:26 $

% Copyright (C) 2005  Darren L. Weber
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

% Created: 08/2005, Darren.Weber_at_radiology.ucsf.edu
%                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ver = '$Revision: 1.1 $ $Date: 2009-01-30 03:49:26 $';
fprintf('\nCTF_CHANNEL_METRICS [v %s]\n',ver(11:15));

meg_sens = ctf_channel_select(ctf,'meg');

coil_inner = zeros(length(meg_sens), 3);
coil_outer = coil_inner;
coilOrient_inner = coil_inner;
coilOrient_outer = coil_inner;

for i = 1:length(meg_sens),
    
    index = meg_sens(i);
    coil_i = ctf.sensor.info(index).hcoil(1).position;
    coil_i = [ coil_i.x; coil_i.y; coil_i.z ];
    coil_o = ctf.sensor.info(index).hcoil(2).position;
    coil_o = [ coil_o.x; coil_o.y; coil_o.z ];
    coil_inner(i,:) = coil_i';
    coil_outer(i,:) = coil_o';
    
    coilOrient_i = ctf.sensor.info(index).hcoil(1).orient;
    coilOrient_i = [ coilOrient_i.x; coilOrient_i.y; coilOrient_i.z ];
    coilOrient_o = ctf.sensor.info(index).hcoil(2).orient;
    coilOrient_o = [ coilOrient_o.x; coilOrient_o.y; coilOrient_o.z ];
    coilOrient_inner(i,:) = coilOrient_i';
    coilOrient_outer(i,:) = coilOrient_o';
    
end

surf_inner.vertices = coil_inner;
surf_inner.faces = convhulln(coil_inner);

surf_outer.vertices = coil_outer;
surf_outer.faces = convhulln(coil_outer);

figure; hold on
patch('faces',surf_inner.faces,'vertices',surf_inner.vertices,...
    'facecolor',[0.9 0.5 0.5],'edgecolor','none',...
    'facealpha',0.8);
patch('faces',surf_outer.faces,'vertices',surf_outer.vertices,...
    'facecolor',[0.5 0.5 0.9],'edgecolor','none',...
    'facealpha',0.2);
legend('inner','outer')
scatter3(coil_inner(:,1), coil_inner(:,2), coil_inner(:,3), 10, 'r', 'filled')
scatter3(coil_outer(:,1), coil_outer(:,2), coil_outer(:,3), 10, 'b', 'filled')

X = coil_outer(:,1);
Y = coil_outer(:,2);
Z = coil_outer(:,3);
U = coilOrient_outer(:,1);
V = coilOrient_outer(:,2);
W = coilOrient_outer(:,3);
quiver3(X, Y, Z, U, V, W, 0)

set(gca,'DataAspectRatio',[1,1,1]);
rotate3d on

coil_separation = coil_outer - coil_inner;
coil_distance = zeros(length(coil_separation),1);
for i = 1:length(coil_separation),
    coil_distance(i) = norm(coil_separation(i,:));
end

stats.coil_distance = coil_distance;
stats.mean = mean(coil_distance);
stats.std  = std(coil_distance);

return
