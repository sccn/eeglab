% showclusts() - Visualize getclusts() clusters. Opens several figures:
%                -  A figure with sorted silhouette values of each cluster
%                -  A figure showing clusters in 3-D, with brightness of
%                     each point proportional to its silhouette value
%                -  Multiple 2-D plots, each for one cluster. Background
%                     color shows interpolated silhouette value. A bright
%                     background means the component fits well in the cluster;
%                     a dark background means it belongs near equally
%                     to another cluster or to no cluster.
% Usage:
%              >>  mutual_info = mi_pairs(activities);
%              >>  clusters = getclusts(mutual_info);
%        Then
%              >>  showclusts(clusters)
%        else  >>  showclusts(clusters,compplotfunc)
% Input:
%
% clusters:    - cluster structure created by getclusts()
%
% Optional Input:
%
% compplotfunc - function call to plot a graphic represenation
%                for each component (for example, 2-D or 3-D
%                scalp maps of ICA components). See example below.
%                The function should use the component number (below, x)
%                as its first input parameter {default: Plot a circle
%                for each component with radius proportional to
%                the component's silhouette value in the cluster}.
% Example:
%
%  % Plot mutual information clusters by plotting component maps
%  %    in the 2-D MDS projection.
%  >>   plotTopo = @(x) topoplot(EEG.icawinv(:,x), EEG.chanlocs,'electrodes','off');
%  >>   showclusts(EEG.miclust,plotTopo);
%
% See also: getclusts(), eeg_miclust(), mi_pairs(), mdscale(), silhouette()
%
% Author: Nima Bigdely Shamlo, SCCN/INC/UCSD, 2006

% Copyright (C) 2006 Nima Bigdely Shamlo, SCCN/INC/UCSD, nima@sccn.ucsd.edu
%
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

function showclusts(miclusters,plotComponentFunction, showAdditionalFigures ,mapSizeCoeff)

if nargin < 3
    showAdditionalFigures = true;
end;

if nargin < 4
    mapSizeCoeff = 1;
end;

if ~isstruct(miclusters)
    error('miclusters input must be a getclusts() output structure');
end

t2d = miclusters.coord2d;
t3d = miclusters.coord3d;
t_manyd = miclusters.coord_manyd;
IDX=miclusters.clusterIDs;
number_of_clusters = miclusters.number_of_clusters;
components = miclusters.allcomponents;

t2d = [t2d components'];        % put component number in the third row
t2d = [t2d IDX];                % put cluster number in the fourth row
t2d = [t2d (1:size(t2d,1))'];   % put original row index number in the fifth row
t3d = [t3d (1:size(t2d,1))'];


%
% Plot silhouette
%

SilhouetteFigureHandle = figure;
[silh h] = silhouette(t_manyd,IDX);
set(gcf,'Name','Cluster Silhouette');
mapAxisSize = 0.09 * mapSizeCoeff;
silhouetteColorSaturation = 3;
silhColor3D= 0.3+((silh*2+1)/2);
silhColor3D = max(0.5,silhColor3D);
silhColor3D = min(1,silhColor3D);
if ~showAdditionalFigures
    close(SilhouetteFigureHandle);
end;

%
% Define cluster colors
%
clusterColormap = colorDeSaturate([0 0 1; 1 0 0; 0 1 0; 1 1 0;1 0 1]);
if number_of_clusters<=5
    cluster_colors=clusterColormap;
else
    cluster_colors=jet(number_of_clusters);
end;

colors = [];
for i=1:size(t2d,1)
    colors = [colors; cluster_colors(IDX(i),:)];
end;

%
%  Plot 3-D Clusters
%
figure3Dhandle = figure('Name','Clusters in 3D','NumberTitle','off');
sphereSize = max(max(t3d(:,1:3))-min(t3d(:,1:3)))/38;

hold on;
legends = '1';
for i=1:number_of_clusters
    t3d_forcluster = t3d(IDX==i,:);
    legends = char(legends,['Cluster ' num2str(i)]);

    %    micluster.cluster{i} =
    %   scatter3(t3d_forcluster(:,1),t3d_forcluster(:,2),t3d_forcluster(:,3),50,IDX(find(IDX==i),:),'filled');
    %     scatter3(t3d_forcluster(:,1),t3d_forcluster(:,2),t3d_forcluster(:,3),50,colors(find(IDX==i),:),'filled');
    %   plot3(t3d_forcluster(:,1),t3d_forcluster(:,2),t3d_forcluster(:,3),'o','MarkerSize',12,'MarkerFaceColor',[.49 1 .63]);

    for j=1:size(t3d_forcluster,1)
        patchHandle(i) =  putSpehere(t3d_forcluster(j,1),...
            t3d_forcluster(j,2),...
            t3d_forcluster(j,3),...
            colors(find(IDX==i,1),:) ...
            * silhColor3D(t3d_forcluster(j,4)),sphereSize);
    end;
%end;

%
% Scale 2-D coordinates to figure dimensions ([0 1)
%
%   t2d_for_axes = t2d(IDX==i,:);
if size(t3d_forcluster,1)>1
    p = pdist(miclusters.cluster{i}.coord2d);
    if isempty(p) || min(p)/max(p)>.05 % check if some are too close. If so,
        % use global 2-D coordinates instead
        t2d_for_axes = [miclusters.cluster{i}.coord2d ...
            miclusters.cluster{i}.components' ...
            repmat(i,size(miclusters.cluster{i}.coord2d,1))];
        t2d_for_axes(:,1) = t2d_for_axes(:,1) - mean(t2d_for_axes(:,1));
        % place their centroid at (0,0)
        t2d_for_axes(:,2) = t2d_for_axes(:,2) - mean(t2d_for_axes(:,2));
    else
        t2d_for_axes = t2d(IDX==i,:);
    end;

    if size(t2d_for_axes,1) > 1
        max_range= max((max(t2d_for_axes(:,1)) - min(t2d_for_axes(:,1))),...
            (max(t2d_for_axes(:,2)) - min(t2d_for_axes(:,2))));
        t2d_for_axes(:,1) = t2d_for_axes(:,1)/ max_range; % set the range in (0,1)
        t2d_for_axes(:,2) = t2d_for_axes(:,2)/ max_range;
        max_range= max((max(t2d_for_axes(:,1)) - min(t2d_for_axes(:,1))),...
            (max(t2d_for_axes(:,2)) - min(t2d_for_axes(:,2))));
    end;

    t2d_for_axes(:,1) = t2d_for_axes(:,1) - min(t2d_for_axes(:,1));
    % place their centroid at (0,0)
    t2d_for_axes(:,2) = t2d_for_axes(:,2) - min(t2d_for_axes(:,2));

    t2d_for_axes(:,1) = 0.5 +  ( t2d_for_axes(:,1) - mean(t2d_for_axes(:,1)));
    % center on (0.5,0.5) of the figure;
    t2d_for_axes(:,2) = 0.5 +  (t2d_for_axes(:,2) - mean(t2d_for_axes(:,2)));

    x = t2d_for_axes(:,1);
    t2d_for_axes(find(x<0.1),1) = 0.1;
    t2d_for_axes(find(x>0.9),1) = 0.9;
    y = t2d_for_axes(:,2);
    t2d_for_axes(find(y<0.1),2) = 0.1;
    t2d_for_axes(find(y>0.9),2) = 0.9;
else
    t2d_for_axes = [.5 .5 miclusters.cluster{i}.components' i i];
end;

positions =t2d_for_axes(:,1:2);
relaxedPositions = relaxCoordinates(positions, mapAxisSize);
t2d_for_axes(:,1:2) = relaxedPositions;

comp_cluster{i} =t2d_for_axes;
miclusters.cluster{i}.t2d_for_axes = comp_cluster{i};


end; %?

hold off;
axis equal;
view(-60,60);
set(gca,'projection','perspective');set(gca,'box','on');
set(gca,'xtick',[]);set(gca,'ytick',[]);set(gca,'ztick',[]);
legends(1,:) = [];

legend(patchHandle,legends);
camlight headlight;
lighting phong;
material shiny;

h = axes('position',[0.42 0.95 0.2 0.09],'Visible','off');
text(.15,.15,['Mutual Information Clusters']);

if ~showAdditionalFigures
    close(figure3Dhandle);
end;

%
% Plot 2D Clusters
%
for cluster_counter = 1:number_of_clusters

    fg = figure('Name',['Cluster ' num2str(cluster_counter)],'NumberTitle','off');
    %
    % adding interpolated silhoette background image
    %
    figureColor = [0.9300    0.9600    1.0000];
    set(gcf,'Color',figureColor);
    h = axes('position',[0 0 1 1],'Visible','off');
    [XI,YI] = meshgrid(0:0.005:1,0:0.005:1);
    if length(find(IDX==cluster_counter))>1
        try
        ZI = griddata(comp_cluster{cluster_counter}(:,1),comp_cluster{cluster_counter}(:,2),...
            silh(find(IDX==cluster_counter)),XI,YI,'v4');
        catch
            ZI = ones(size(XI,1),size(XI,2)) * silh(find(IDX==cluster_counter, 1));
        end;
    else
        ZI = ones(size(XI,1),size(XI,2)) * silh(find(IDX==cluster_counter));
    end;
    axis tight;axis ij;
    z= (ZI*silhouetteColorSaturation+1)/2;  z = max(0,z); z = min(1,z);
    im=zeros(size(z,1),size(z,2),3);
    im(:,:,1) = z*figureColor(1);im(:,:,2) = z*figureColor(2);im(:,:,3) = z*figureColor(3);
    im = flipdim(im,1);
    image([min(XI(:)),max(XI(:))],[min(YI(:)),max(YI(:))],im);  axis off;axis tight;

    for i=1:size(comp_cluster{cluster_counter},1)
        % change axis size according to silhuoette value
        silhAxisSize= (silh(comp_cluster{cluster_counter}(i,5))*.5+1)/2;
        silhAxisSize = max(0.3,silhAxisSize);
        silhAxisSize = min(1,silhAxisSize);

        % for now only, an additional parajeter should be assigned for this
        % later
        if showAdditionalFigures
            mapAxisSize = 0.09 * silhAxisSize*1.5 * mapSizeCoeff;
        end;

        rect = [comp_cluster{cluster_counter}(i,1) - mapAxisSize/2 ...
            comp_cluster{cluster_counter}(i,2) - mapAxisSize/2 mapAxisSize mapAxisSize];
        h = axes('position',rect);
        axis tight;
        silhColor= 1-((silh(comp_cluster{cluster_counter}(i,5))*silhouetteColorSaturation+1)/2);
        silhColor = max(0,silhColor); silhColor = min(1,silhColor);
        % silhColor = 1 - ...
        %             im(round((-mapAxisSize/2 + comp_cluster{cluster_counter}(i,1))*size(im,1)),...
        %             round((comp_cluster{cluster_counter}(i,2))*size(im,2)));
        if silhColor <.5
            silhColor= .3;
        else
            silhColor = 1;
        end;

        % if showHead
        %   headplot(EEG.icawinv(:,comp_cluster{cluster_counter}(i,3)),spline_file ,...
        %                'electrodes','off');
        % else
        %   tmpobj = topoplot(EEG.icawinv(:,comp_cluster{cluster_counter}(i,3)), ...
        %               EEG.chanlocs);
        % end;

        figure(fg);
        if nargin>1
            if nargin(plotComponentFunction)==1
                plotComponentFunction(comp_cluster{cluster_counter}(i,3));
            else
                plotComponentFunction(comp_cluster{cluster_counter}(i,3),cluster_counter);
            end;
            title(num2str(comp_cluster{cluster_counter}(i,3)),'color',silhColor*[1 1 1]);
            %,'BackgroundColor',(1-silhColor)*[1 1 1]
        else
            set(h,'visible','off');
            rectangle('Position',[0,0,1,1],'Curvature',[1,1],...
                'FaceColor',cluster_colors(cluster_counter,:));

            % text(.5,.4,num2str(comp_cluster{cluster_counter}(i,3)),...
            %        'FontSize',14,'color',imcomplement(cluster_colors(cluster_counter,:)),...
            %        'HorizontalAlignment','center','VerticalAlignment','middle');

            annotation('textbox',[comp_cluster{cluster_counter}(i,1) - mapAxisSize/2 ...
                comp_cluster{cluster_counter}(i,2) - mapAxisSize/2 ...
                mapAxisSize mapAxisSize],...
                'string',num2str(comp_cluster{cluster_counter}(i,3)),...
                'color',imcomplement(cluster_colors(cluster_counter,:)),...
                'LineStyle','none',...
                'HorizontalAlignment','center',...
                'VerticalAlignment','middle');
            daspect([1,1,1])
        end;

        axes(h);
        drawnow;
    end;

    h = axes('position',[0.4 0.9 0.2 0.09],'Visible','off');
    text(.15,.15,...
        ['Cluster ' num2str(cluster_counter)],...
        'BackgroundColor',cluster_colors(cluster_counter,:),...
        'color',imcomplement(cluster_colors(cluster_counter,:)));
    axcopy;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Sub-functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hs = putSpehere(x,y,z,c,sphereSize)
[Xs Ys Zs]=sphere(20);
hs = surf((Xs*sphereSize)+x, ...
    (Ys*sphereSize)+y, ...
    (Zs*sphereSize)+z);
set(hs,'EdgeColor','none', ...
    'FaceColor',c, ...
    'FaceLighting','phong', ...
    'AmbientStrength',0.4, ...
    'DiffuseStrength',0.8, ...
    'SpecularStrength',0.1, ...
    'SpecularExponent',20, ...
    'BackFaceLighting','lit');

function desatcolor = colorDeSaturate(color)
maxSaturation = .85;
hsv = rgb2hsv(color);
t = hsv(:,2);
t(t>maxSaturation) =maxSaturation;
hsv(:,2) = t;
desatcolor = hsv2rgb(hsv);

% push overlapping ICs apart
function coord = relaxCoordinates(coord, radius)
% if there is more than one point
if size(coord, 1) > 1
    counter = 1;
    while counter ==1 | (length(find((displacementLengths(:)>1e-9))>0) && counter < 5000)
        distances = squareform(pdist(coord));
        displacementLengths = 0.01 * (radius - distances) / radius;
        displacementLengths(displacementLengths<0) = 0;
        displacementLengths = displacementLengths - diag(diag(displacementLengths));
        
        xs = repmat(coord(:,1),1,size(coord,1));
        ys = repmat(coord(:,2),1,size(coord,1));
        xDifferences = xs - xs';
        yDifferences = ys - ys';
        
        xChanges = xDifferences .* displacementLengths;
        yChanges = yDifferences .* displacementLengths;
        
        
        for i=1:size(coord,1)
            coord(i,1) = coord(i,1) + sum(xChanges(i,:));
            coord(i,2) = coord(i,2) + sum(yChanges(i,:));
        end;
        counter = counter + 1;
    end;
end;
