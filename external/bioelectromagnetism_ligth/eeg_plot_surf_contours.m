function [ p ] = eeg_plot_surf_contours(p,mode)

% eeg_plot_surf_contours - plot 3D contours on triangulated surface
% 
% Usage: [p] = eeg_contours_3d(p,[mode])
% 
% p is the eeg_toolbox struct (see eeg_toolbox_defaults)
%
% mode is:
%
% 'rb' for +ve red, -ve blue
% 'bw' for +ve solid black, -ve dashed black (default)
% 
% In development!
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:52 $

% Licence:  Gnu GPL, no express or implied warranties
% History:  04/03 Darren.Weber_at_radiology.ucsf.edu
%           obtained permission from Robert Oostenveld to
%           distribute his code below under the GPL
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('\nEEG_PLOT_SURF_CONTOURS (%s)...\n',['$Revision: 1.1 $']); tic;

% confirm that a patch surface is available
if isfield(p,'mesh'),
    if isfield(p.mesh,'data'),
        if isfield(p.mesh.data,'timeseries'),
            if isempty(p.mesh.data.Cdata{p.mesh.current}),
                msg = sprintf('...p.mesh.data.Cdata{%d} is empty\n',p.mesh.current);
                error(msg);
            end
        end
    end
else
    error('...p.mesh.data is empty - load mesh first\n');
end

if ~exist('mode','var'), mode = 'bw'; end
if isempty(mode),        mode = 'bw'; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine/validate the min, max surface color data range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select mesh timeseries values at samplePoint
switch p.mesh.data.meshtype{p.mesh.current},
case {'scalp','elec'},
    samplePoint = p.volt.samplePoint;
otherwise
    samplePoint = p.mesh.samplePoint;
end

% get vertices and faces of mesh
vertices = p.mesh.data.vertices{p.mesh.current};
faces = p.mesh.data.faces{p.mesh.current};

% get number of vertices
nvert = size(p.mesh.data.vertices{p.mesh.current},1);
% Assume more vertices than time points
[s1,s2] = size(p.mesh.data.Cdata{p.mesh.current});
if isequal(nvert,s1), % vertices in rows, timepoints in columns
    meshCdata = p.mesh.data.Cdata{p.mesh.current}(:,samplePoint);
else,                 % vertices in columns, timeseries in rows
    meshCdata = p.mesh.data.Cdata{p.mesh.current}(samplePoint,:)';
end

switch p.rangeMethod
case 'minmaxall', % Min/Max,all points
    fprintf('...estimating color data range, min/max all time points.\n');
    p.maximumIntensity = max(max(p.mesh.data.Cdata{p.mesh.current}));
    p.minimumIntensity = min(min(p.mesh.data.Cdata{p.mesh.current}));
case 'minmaxone', % Min/Max, single point
    fprintf('...estimating color data range, min/max single time point.\n');
    % get number of vertices
    p.maximumIntensity = max(max(meshCdata));
    p.minimumIntensity = min(min(meshCdata));
case 'minmaxabs', % Min/Max, Absolute
    fprintf('...estimating color data range, abs min/max single time point.\n');
    absmax = max(max(abs(meshCdata)));
    p.maximumIntensity =  absmax;
    p.minimumIntensity = -absmax;
otherwise
    % check that specified intensity range is defined
    fprintf('...checking predefined color data range.\n');
    if isempty(p.maximumIntensity),
        fprintf('...estimating color data range, min/max single time point.\n');
        p.maximumIntensity = max(max(meshCdata)); end
    if isempty(p.minimumIntensity),
        fprintf('...estimating color data range, min/max single time point.\n');
        p.minimumIntensity = min(min(meshCdata)); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the contour values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch p.contour.stepMethod
case 0  % StepSize method
    fprintf('...using contour step size: %6.2f\n',p.contour.stepSize);
    p.contour.levels = eeg_contour_levels( p.contour.stepSize, meshCdata );
    p.contour.Nsteps = length(p.contour.levels);
otherwise % Number of steps method
    fprintf('...using number of contours: %d\n',p.contour.Nsteps);
    p.contour.stepSize = abs(p.maximumIntensity - p.minimumIntensity) / p.contour.Nsteps;
    p.contour.levels = eeg_contour_levels( p.contour.stepSize, meshCdata );
end


% Exit if there are no contour levels
if isempty(p.contour.levels),
    fprintf('...no contours in this data range.\n');
    % remove any current contours
    if isfield(p.contour,'patches'),
        if ~isempty(p.contour.patches),
            handleIndex = find(ishandle(p.contour.patches));
            delete(p.contour.patches(handleIndex));
        end
    end
    p.contour.patches = [];
    return
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code below adapted from triplot.m (c) 2001 Robert Oostenveld ...
% This code here is licenced under the GPL, as of April 8th 2003.

fprintf('...calculating surface contours.\n');


% create an Nx3 matrix of color data values at each face vertex
FaceCdata = meshCdata(faces);
% find the min/max vertex values of each face
FaceCdata_min = min(FaceCdata, [], 2);
FaceCdata_max = max(FaceCdata, [], 2);

for contourIndex=1:length(p.contour.levels),
    
    % Find all the faces containing any vertices within the contour level
    contourLevel = p.contour.levels(contourIndex);
    use = contourLevel>=FaceCdata_min & contourLevel<=FaceCdata_max;
    faceIndices = find(use)';
    
    intersect1 = [];
    intersect2 = [];
    
    counter = 0;
    for faceIndex=faceIndices,
        
        xyz  = vertices(faces(faceIndex,:), :); % 3x3, rows = vertexIJK, columns = XYZ
        v = FaceCdata(faceIndex,:);             % 1x3, Cdata at rows of xyz
        
        % find line through face that passes through contourLevel ?
        
        la(1) = (contourLevel-v(1)) / (v(2)-v(1));	% abscissa between vertex 1 and 2
        la(2) = (contourLevel-v(2)) / (v(3)-v(2));	% abscissa between vertex 2 and 3
        la(3) = (contourLevel-v(3)) / (v(1)-v(3));	% abscissa between vertex 1 and 2
        
        abc(1,:) = xyz(1,:) + la(1) * (xyz(2,:) - xyz(1,:));
        abc(2,:) = xyz(2,:) + la(2) * (xyz(3,:) - xyz(2,:));
        abc(3,:) = xyz(3,:) + la(3) * (xyz(1,:) - xyz(3,:));
        
        counter = counter + 1;
        select  = find(la>=0 & la<=1);
        
        intersect1(counter, :) = abc(select(1),:);
        intersect2(counter, :) = abc(select(2),:);
    end
    
    % remember the details for external reference
    contour(contourIndex).level = contourLevel;
    contour(contourIndex).n     = counter;
    contour(contourIndex).intersect1 = intersect1;
    contour(contourIndex).intersect2 = intersect2;
end


% collect all different contourlevels
intersect1 = [];
intersect2 = [];
contourLevel = [];
for contourIndex=1:length(p.contour.levels)
    intersect1 = [intersect1; contour(contourIndex).intersect1];
    intersect2 = [intersect2; contour(contourIndex).intersect2];
    contourLevel = [contourLevel; ones(contour(contourIndex).n,1) * p.contour.levels(contourIndex)];
end

X = [intersect1(:,1) intersect2(:,1)]';
Y = [intersect1(:,2) intersect2(:,2)]';
C = [contourLevel(:) contourLevel(:)]';

if size(vertices,2)>2
    % 3D contours
    Z = [intersect1(:,3) intersect2(:,3)]';
else
    % 2D contours
    Z = zeros(2, length(contourLevel));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the surface contours
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove any current contours
if isfield(p.contour,'patches'),
    if ~isempty(p.contour.patches),
        handleIndex = find(ishandle(p.contour.patches));
        delete(p.contour.patches(handleIndex));
    end
end

p.contour.patches = [];

fprintf('...plotting surface contours.\n');

for i=1:length(contourLevel)
    
    % make colour contours (default)
    h = patch('Tag','CONTOUR','XData', X(:,i), 'Ydata', Y(:,i), ...
        'ZData', Z(:,i), 'CData', C(:,i), ...
        'facecolor','none','edgecolor','flat', ...
        'linestyle', '-', 'linewidth', 1, ...
        'userdata',contourLevel(i));
    
    switch mode,
    case 'rb',
        % make red-blue contours
        if     contourLevel(i)>0,   edgecolor = 'red';
        elseif contourLevel(i)<0,   edgecolor = 'blue';
        else                        edgecolor = 'black';
        end
    otherwise
        % make black-white contours
        if     contourLevel(i)>0,   edgecolor = [.4 .4 .4];
        elseif contourLevel(i)<0,   edgecolor = [.2 .2 .2];
        else                        edgecolor = [.5 .5 .5];
        end
    end
    % apply contour styles
    if     contourLevel(i)>0,   set(h,'edgecolor',edgecolor,'linestyle', '-');
    elseif contourLevel(i)<0,   set(h,'edgecolor',edgecolor,'linestyle', ':','linewidth', 1);
    else                        set(h,'edgecolor',edgecolor,'linewidth', 1.5);
    end
    
    p.contour.patches = [p.contour.patches; h];
end

t = toc;
fprintf('...done (%6.2f sec).\n\n',t);

return
