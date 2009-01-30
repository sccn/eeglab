function [Handles] = ctf_channel_plot(ctf,surf,CHAN,Handles)

% ctf_channel_plot - plot ctf.sensor data
%
% [Handles] = ctf_channel_plot(ctf [,surf][,CHAN][,H])
%
% ctf  - a struct returned by ctf_read
% surf - a surface tesselation struct, surf.vertices and surf.faces;
%        the default is to create a convex hull of the sensor locations
% CHAN - see ctf_channel_select for options (default, CHAN = 'meg')
%
% Handles - handles struct, eg:
%           H.Hf - figure handle to plot into (default, H.Hf = figure;)
%           H.Hs - sensor handles
%           H.Hp - patch handles
%           H.Hl - light handles
%
% All inputs, except ctf, are optional.  An empty input [] will be replaced
% with the default.
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

% Copyright (C) 2004  Darren L. Weber
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

% Modified: 05/2004, Darren.Weber_at_radiology.ucsf.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('ctf','var'), error('no ctf input'); end
if isempty(ctf), error('no ctf input'); end

% check the figure handle for initialization
if ~exist('Handles','var'), 
  Handles.Hf = figure;
  colordef(Handles.Hf,'white')
end
if isempty(Handles), 
  Handles.Hf = figure;
  colordef(Handles.Hf,'white')
end
if ~isfield(Handles,'Hf'),
  Handles.Hf = figure;
  colordef(Handles.Hf,'white')
end
if isempty(Handles.Hf),
  Handles.Hf = figure;
  colordef(Handles.Hf,'white')
end

% check other fields of the Handles
if ~isfield(Handles,'Hp'), Handles.Hp = {}; end
if ~isfield(Handles,'Hs'), Handles.Hs = {}; end
if ~isfield(Handles,'Hl'), Handles.Hl = []; end

usesurf = 1; % assume input surf
if ~exist('surf','var'), usesurf = 0; end
if isempty(surf), usesurf = 0; end

if ~exist('CHAN','var'), CHAN = 'meg'; end
if CHAN,
  % Get the channel indices to plot
  [CHAN,type] = ctf_channel_select(ctf,CHAN);
  sensors.vertices = zeros(length(CHAN),3);
  for i = 1:length(CHAN),
    sensors.vertices(i,:) = ctf.sensor.info(CHAN(i)).location(:,1)';
  end
end

ver = '$Revision: 1.1 $ $Date: 2009-01-30 03:49:26 $';
fprintf('\nCTF_CHANNEL_PLOT [v %s]\n',ver(11:15));

figure(Handles.Hf); hold on
set(Handles.Hf,'Renderer','OpenGL')

if CHAN,
  % plot the sensor cloud
  Handles.Hs{end+1} = scatter3(sensors.vertices(:,1),...
                               sensors.vertices(:,2),...
                               sensors.vertices(:,3),...
                               50,'bo','filled');
end

if ~usesurf,
  % nothing provided, assume we want
  % to create a convex hull of the sensors
  surf.vertices = sensors.vertices;
  surf.faces = convhulln(surf.vertices);
end
Nvert = length(surf.vertices);

% assume surface is a skin color
skinRGB = [1, 0.85, 0.75];
vertexColors = ones(Nvert,1) * skinRGB;

% plot the surface patch (eg, sensors or scalp)
Handles.Hp{end+1} = patch('faces',surf.faces,'vertices',surf.vertices,...
                          'facecolor','interp',...
                          'FaceVertexCData', vertexColors,...
                          'edgecolor', 'none',...
                          'facealpha',0.5);

% view parameters
view([40,20])
daspect([1 1 1])
axis equal
axis tight
grid on

lighting phong
material dull

% only add a light if it doesn't exist already
if ~isfield(Handles,'Hl'),
  Handles.Hl = camlight;
elseif isempty(Handles.Hl),
  Handles.Hl = camlight;
end

% The CTF MEG Head Coordinate System is
% +X anterior, +Y left, +Z superior
xlabel('x (cm)')
ylabel('y (cm)')
zlabel('z (cm)')

rotate3d on

return
