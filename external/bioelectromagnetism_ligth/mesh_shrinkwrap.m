function [FV, Edges] = mesh_shrinkwrap(vol,FV,smooth,vthresh,interpVal,...
    fitval,fittol,fititer,fitchange,fitvattr)

% mesh_shrinkwrap - Tesselate the surface of a 3D volume
% 
% [FV, Edges] = mesh_shrinkwrap(vol,FV,smooth,vthresh,interpVal,...
%               fitval,fittol,fititer,fitchange,fitvattr)
% 
% vol       - a 3D image volume
% FV        - input tesselation; if empty, sphere tesselation 
%             is created.  FV has fields FV.vertices, FV.faces
% smooth    - Gaussian smoothing of vol (5x5x5 kernel), 0|1 (default 0)
% vthresh   - binarise vol at fitval, 0|1 (default 0)
% interpVal - use radial interpolation (faster) or incremental
%             radial shrink (slower), 0|1 (default 1, faster)
% fitval    - image intensity to shrink wrap (default 20)
% fittol    - image intensity tolerance (default 5)
% fititer   - max number of iterations to fit (default 200)
% fitchange - least sig. change in intensity values 
%             between fit iterations (default 2)
% fitvattr  - vertex attraction (constraint), 0:1, smaller
%             values are less constraint; close to 0 for 
%             no constraint is useful when dealing with 
%             binary volumes, otherwise 0.4 (40%) seems OK
% 
% FV        - a struct with 2562 vertices and 5120 faces
% Edges     - a [2562,2562] matrix of edge connectivity for FV
% 
% An alternative to isosurface for volumes with an external 
% surface that can be shrink-wrapped. It has been developed to
% find the scalp surface for MRI of the human head.
% 
% It starts with a sphere tesselation (large radius) and moves
% each vertex point toward the centre of the volume until it
% lies at or near the fitval.  The function is not optimised
% for speed, but it should produce reasonable results.
% 
% Example of creating a scalp tesselation for SPM T1 MRI template:
% 
%   avw = avw_read('T1');
%   FV = mesh_shrinkwrap(avw.img,[],0,0,intensity,5.0,50,0.5,0.4);
%   patch('vertices',FV.vertices,'faces',FV.faces,'facecolor',[.6 .6 .6]);
% 
% Example of creating a skull from FSL BET skull volume:
%   
%   avw = avw_read('T1_skull');
%   FV = mesh_shrinkwrap(avw.img,[],1,0,intensity,0.2,10,0.005,0.1);
%   patch('vertices',FV.vertices,'faces',FV.faces,'facecolor',[.6 .6 .6]);
% 
% ***** NOTE *****
% An important limitation at present is that it doesn't
% read the image dimensions, but assumes they are 1mm^3.  Further
% versions might take a standard analyze volume and check the
% header to scale the result according to the image dimensions. At
% present, this must be done with the results of this function.  The
% output vertex coordinates are in mm with an origin at (0,0,0), 
% which lies at the center of the MRI volume.
% 
% See also: ISOSURFACE, SPHERE_TRI, MESH_REFINE_TRI4,
%           MESH_BEM_SHELLS_FUNC, MESH_BEM_SHELLS_SCRIPT
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:57 $

% Licence:  GNU GPL, no implied or express warranties
% History:  07/2002, Darren.Weber_at_radiology.ucsf.edu
%                    - created to provide alternative to
%                      isosurface in matlab R12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

version = '[$Revision: 1.1 $]';
fprintf('MESH_SHRINKWRAP [v%s]\n',version(12:16));  tic;

% Parse arguments

if ~exist('vol','var'), error('...no input volume\n');
elseif isempty(vol),    error('...empty input volume\n');
end
if isstruct(vol),
    if isfield(vol,'img'), vol = vol.img;
    else
        error('...input volume is a struct, it must be a 3D image volume only\n');
    end
end

if ~exist('smooth','var'),   smooth = 0;
elseif isempty(smooth),      smooth = 0;
end

if ~exist('vthresh','var'),  vthresh = 0;
elseif isempty(vthresh),     vthresh = 0;
end

if ~exist('interpVal','var'), interpVal = 1;
elseif isempty(interpVal),    interpVal = 1;
end

if ~exist('fitval','var'),   fit.val = 20;
elseif isempty(fitval),      fit.val = 20;
else                         fit.val = fitval;
end

if ~exist('fittol','var'),   fit.tol = 5;
elseif isempty(fittol),      fit.tol = 5;
else                         fit.tol = fittol;
end

if fit.val <= fit.tol,
    error('...must use fit tolerance < fit value\n');
end

if ~exist('fititer','var'),  fit.iter = 200;
elseif isempty(fititer),     fit.iter = 200;
else                         fit.iter = fititer;
end

if ~exist('fitchange','var'),fit.change = 2;
elseif isempty(fitchange),   fit.change = 2;
else                         fit.change = fitchange;
end

if ~exist('fitvattr','var'), fit.vattr = 0.4;
elseif isempty(fitvattr),    fit.vattr = 0.4;
else                         fit.vattr = fitvattr;
end
if fit.vattr > 1,
    fprintf('...fit vertattr (v) must be 0 <= v <= 1, setting v = 1\n');
    fit.vattr = 1;
end
if fit.vattr < 0,
    fprintf('...fit vertattr (v) must be 0 <= v <= 1, setting v = 0.\n');
    fit.vattr = 0;
end

% MAIN

% Find approximate centre of volume
xdim = size(vol,1);
ydim = size(vol,2);
zdim = size(vol,3);

origin(1) = floor(xdim/2);
origin(2) = floor(ydim/2);
origin(3) = floor(zdim/2);

% Check whether to create a sphere tesselation
% or use an input tesselation as the start point
sphere = 0;
if ~exist('FV','var'),
    sphere = 1;
elseif ~isfield(FV,'vertices'),
    sphere = 1;
elseif ~isfield(FV,'faces'),
    sphere = 1;
elseif isempty(FV.vertices),
    sphere = 1;
elseif isempty(FV.faces),
    sphere = 1;
end
if sphere,
    % Create a sphere tesselation to encompass the volume
    radius = max([xdim ydim zdim]) / 1.5;
    FV = sphere_tri('ico',4,radius); % 2562 vertices
else
    fprintf('...using input FV tesselation...\n');
end


% the 'edge' matrix is the connectivity of all vertices,
% used to find neighbours during movement of vertices,
% including smoothing the tesselation
FV.edge = mesh_edges(FV);


% Shift the centre of the sphere to the centre of the volume
centre = repmat(origin,size(FV.vertices,1),1);
FV.vertices = FV.vertices + centre;


% 'Binarise' the volume, removing all values below
% a threshold, setting all others to threshold
if vthresh,
    fprintf('...thresholding volume...'); tic;
    Vindex = find(vol < fit.val);
    vol(Vindex) = 0;
    Vindex = find(vol >= fit.val);
    vol(Vindex) = fit.val;
    t = toc; fprintf('done (%5.2f sec)\n',t);
end
binvol = find(vol > 1);

% smooth the volume
if smooth,
	fprintf('...gaussian smoothing (5-10 minutes)...'); tic;
	vol = smooth3(vol,'gaussian',5,.8);
	t = toc; fprintf('done (%5.2f sec)\n',t);
end


% Now begin recursion
fprintf('...fitting...\n');	tic;

i = 1;
Fminima = 0;
intensityAtDMean = [0 0];

while i <= fit.iter,
    
    if interpVal,
        % use radial interpolation method, moving directly
        % to the intensity value nearest correct intensity
        [FV, intensityAtD, D] = locate_val(FV,vol,origin,fit);
    else
        % use incremental method, moving along radial line
        % gradually until finding correct intensity
        [FV, intensityAtD, D] = shrink_wrap(FV,vol,origin,fit);
    end
    
    intensityAtDMean(1) = intensityAtDMean(2);
    intensityAtDMean(2) = mean(intensityAtD);
    
    fprintf('...distance:  mean = %8.4f mm, std = %8.4f mm\n',mean(D),std(D));
    fprintf('...intensity: mean = %8.4f,    std = %8.4f\n',...
        mean(intensityAtD),std(intensityAtD));
    fprintf('...real iteration: %3d\n',i);
    
    % Is the mean distance reasonable?
    if mean(D) < 0.5,
        error('...mean distance < 0.5 mm!\n');
    end
    
    % MDifVal is the mean of the absolute difference
    % between the vertex intensity and the fit intensity
    MDifVal = abs(intensityAtDMean(2) - fit.val);
    
    % Is the mean difference within the tolerance range?
    if MDifVal < fit.tol,
        fprintf('...mean intensity difference < tolerance (%5.2f +/- %5.2f)\n',...
            fit.val,fit.tol);
        break;
    else
        fprintf('...mean intensity difference > tolerance (%5.2f +/- %5.2f)\n',...
            fit.val,fit.tol);
    end
    
    % How much has the intensity values changed?
    if (i > 1) & intensityAtDMean(2) > 0,
        if intensityAtDMean(2) - intensityAtDMean(1) < fit.change,
            fprintf('...no significant intensity change (< %5.2f) in this iteration\n',...
                fit.change);
            Fminima = Fminima + 1;
            if Fminima >= 5,
                fprintf('...no significant intensity change in last 5 iterations\n');
                break;
            end
        else
            Fminima = 0;
        end
    end
    
    % Ensure that iterations begin when MDifVal is truly initialised
    if isnan(MDifVal),
        i = 1;
    else,
        i = i + 1;
    end
    
end

FV = mesh_smooth(FV,origin,fit.vattr);

% Remove large edges matrix from FV
Edges = FV.edge;
FV = struct('vertices',FV.vertices,'faces',FV.faces);

% Now center the output vertices at 0,0,0 by subtracting
% the volume centroid
FV.vertices(:,1) = FV.vertices(:,1) - origin(1);
FV.vertices(:,2) = FV.vertices(:,2) - origin(2);
FV.vertices(:,3) = FV.vertices(:,3) - origin(3);

t=toc; fprintf('...done (%5.2f sec).\n\n',t);

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FV, intensityAtD, D] = locate_val(FV,vol,origin,fit),
    
    xo = origin(1); yo = origin(2); zo = origin(3);
    
    Nvert = size(FV.vertices,1);
    progress = round(.1 * Nvert);
    
    % Initialise difference intensity & distance arrays
    intensityAtD = zeros(Nvert,1);
    D    = intensityAtD;
    
    % Find distance and direction cosines for line from 
    % origin to all vertices
    XV = FV.vertices(:,1);
    YV = FV.vertices(:,2);
    ZV = FV.vertices(:,3);
    DV = sqrt( (XV-xo).^2 + (YV-yo).^2 + (ZV-zo).^2 );
    LV = (XV-xo)./DV; % cos alpha
    MV = (YV-yo)./DV; % cos beta
    NV = (ZV-zo)./DV; % cos gamma
    
    % Check for binary volume data, if empty, binary
    binvol = find(vol > 1);
    
    % Locate each vertex at a given fit value
    tic
    for v = 1:Nvert,
        
        if v > progress,
            fprintf('...interp3 processed %4d of %4d vertices',progress,Nvert);
            t = toc; fprintf(' (%5.2f sec)\n',t);
            progress = progress + progress;
        end
        
        % Find direction cosines for line from origin to vertex
        x = XV(v);
        y = YV(v);
        z = ZV(v);
        d = DV(v);
        l = LV(v); % cos alpha
        m = MV(v); % cos beta
        n = NV(v); % cos gamma
        
        % find discrete points between origin
        % and vertex + 50% of vertex distance
        points = 250;
        
        Darray = linspace(0,(d + .2 * d),points);
        
        L = repmat(l,1,points);
        M = repmat(m,1,points);
        N = repmat(n,1,points);
        
        XI = (L .* Darray) + xo;
        YI = (M .* Darray) + yo;
        ZI = (N .* Darray) + zo;
        
        % interpolate volume values at these points
        VI = interp3(vol,YI,XI,ZI,'*linear');
        
        % do we have a binary volume (no values > 1)
        if isempty(binvol),
            maxindex = max(find(VI>0));
            if maxindex,
                D(v) = Darray(maxindex);
            end
        else
            % find the finite values of VI
            index = max(find(VI(isfinite(VI))));
            if index,
                
                % Find nearest volume value to the required fit value
                nearest = max(find(VI >= fit.val));
                
                %[ nearest, value ] = NearestArrayPoint( VI, fit.val );
                
                % Check this nearest index against a differential
                % negative peak value
                %diffVI = diff(VI);
                %if max(VI) > 1,
                %    diffindex = find(diffVI < -20);
                %else
                % probably a binary volume
                %    diffindex = find(diffVI < 0);
                %end
                
                % now set d
                if nearest,
                    D(v) = Darray(nearest);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constrain relocation by fit.vattr, 
        % some % of distance from neighbours
        
        vi = find(FV.edge(v,:));  % the neighbours' indices
        X = FV.vertices(vi,1);    % the neighbours' vertices
        Y = FV.vertices(vi,2);
        Z = FV.vertices(vi,3);
        
        % Find neighbour distances
        DN = sqrt( (X-xo).^2 + (Y-yo).^2 + (Z-zo).^2 );
        % Find mean distance of neighbours
        DNmean = mean(DN);
        
        minattr = fit.vattr;
        maxattr = 1 + fit.vattr;
        
        if D(v) < (minattr * DNmean),
            D(v) = minattr * DNmean;
        end
        if D(v) > (maxattr * DNmean),
            D(v) = maxattr * DNmean;
        end
        if D(v) == 0, D(v) = DNmean; end
        
        % relocate vertex to new distance
        x = (l * D(v)) + xo;
        y = (m * D(v)) + yo;
        z = (n * D(v)) + zo;
        
        FV.vertices(v,:) = [ x y z ];
        
        % Find intensity value at this distance
        intensityAtD(v) = interp1(Darray,VI,D(v),'linear');
        
    end
    
    if isempty(binvol),
        % check outliers and smooth twice for binary volumes
        FV = vertex_outliers(FV, D, origin);
        FV = mesh_smooth(FV,origin,fit.vattr);
    end
    FV = mesh_smooth(FV,origin,fit.vattr);
    
return







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FV, intensityAtD, D] = shrink_wrap(FV,vol,origin,fit),
    
    xo = origin(1); yo = origin(2); zo = origin(3);
    
    Nvert = size(FV.vertices,1);
    
    intensityAtD = zeros(Nvert,1);
    D    = intensityAtD;
    
    for v = 1:Nvert,
        
        x = FV.vertices(v,1);
        y = FV.vertices(v,2);
        z = FV.vertices(v,3);
        
        % Find distance of vertex from origin
        d = sqrt( (x-xo)^2 + (y-yo)^2 + (z-zo)^2 );
        
        % check whether vertex is already at fit.val in vol
        
        volval = vol_val(vol,x,y,z);
        
        if isnan(volval) | (volval < fit.val - fit.tol) | (volval > fit.val + fit.tol),
            
            % Find direction cosines for line from centre to vertex
            l = (x-xo)/d; % cos alpha
            m = (y-yo)/d; % cos beta
            n = (z-zo)/d; % cos gamma
            
            % now modify d by fit.dist
            if isnan(volval) | (volval < fit.val - fit.tol),
                d = d - (d * fit.vattr);
            else
                d = d + (d * fit.vattr);
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Constrain relocation by fit.vattr,
            % some % of distance from neighbours
            
            vi = find(FV.edge(v,:));  % the neighbours' indices
            X = FV.vertices(vi,1);    % the neighbours' vertices
            Y = FV.vertices(vi,2);
            Z = FV.vertices(vi,3);
            
            % Find neighbour distances
            DN = sqrt( (X-xo).^2 + (Y-yo).^2 + (Z-zo).^2 );
            % Find mean distance of neighbours
            DNmean = mean(DN);
            
            minattr = fit.vattr;
            maxattr = 1 + fit.vattr;
            
            if d < (minattr * DNmean),
                d = minattr * DNmean;
            elseif d > (maxattr * DNmean),
                d = maxattr * DNmean;
            end
            
            % locate vertex at this new distance
            x = (l * d) + xo;
            y = (m * d) + yo;
            z = (n * d) + zo;
            
            FV.vertices(v,:) = [ x y z ];
        end
        
        intensityAtD(v) = vol_val(vol,x,y,z);
        D(v) = d;
    end
    
    FV = mesh_smooth(FV,origin,fit.vattr);
    
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Currently not calling this function (Oct 02)
function [val] = vol_val(vol,x,y,z),
    
    % This function just ensures that xyz are
    % actually within the volume before trying
    % to get a volume value
    
    val = nan; % assume zero value
    
    x = round(x);
    y = round(y);
    z = round(z);
    
    if x > 0 & y > 0 & z > 0,
        
        % get bounds of vol
        Xv = size(vol,1);
        Yv = size(vol,2);
        Zv = size(vol,3);
        
        if x <= Xv & y <= Yv & z <= Zv,
            % OK return volume value at xyz
            val = vol(x,y,z);
        end
    end
    
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FV] = vertex_outliers(FV, D, origin),
    
    xo = origin(1); yo = origin(2); zo = origin(3);
    
    % Screen FV for outlying vertices, using
    % mean +/- 2 * stdev of distance from origin
	DistMean = mean(D);
	DistStDev = std(D);
	DistMax = DistMean + 2 * DistStDev;
	DistMin = DistMean - 2 * DistStDev;
	
	for v = 1:size(FV.vertices,1),
        
        if D(v) >= DistMax,
            D(v) = DistMean;
            relocate = 1;
        elseif D(v) <= DistMin,
            D(v) = DistMean;
            relocate = 1;
        else
            relocate = 0;
        end
        
        if relocate,
            x = FV.vertices(v,1);
            y = FV.vertices(v,2);
            z = FV.vertices(v,3);
            
            % Find direction cosines for line from centre to vertex
            d = sqrt( (x-xo)^2 + (y-yo)^2 + (z-zo)^2 );
            l = (x-xo)/d; % cos alpha
            m = (y-yo)/d; % cos beta
            n = (z-zo)/d; % cos gamma
            
            % relocate vertex to this new distance
            x = (l * D(v)) + xo;
            y = (m * D(v)) + yo;
            z = (n * D(v)) + zo;
            
            FV.vertices(v,:) = [ x y z ];
        end
	end
return
