% dipoledensity() - compute and optionally plot a measure of the 3-D spatial 
%                   (in)homogeneity of a specified (large) set of 1- or 2-dipole 
%                   component models, either as physical dipole density or as 
%                   dipole-position entropy across subjects. In either case, 
%                   take into account either all the dipoles, or only the nearest 
%                   dipole from each of the subjects. If no output arguments, 
%                   or if 'plot','on', paints a 3-D density|entropy brain image 
%                   on slices of the Montreal Neurological Institute (MNI) mean 
%                   MR brain image ('standard_BESA/avg152t1.mat'). Calls
%                   dipplot(), 
%                   mri3dplot(), and Fieldtrip function ft_inside_headmodel(). 
% Usage:
%               >> [dens3d mri] = dipoledensity( dipoles, 'key',val, ... );
%
% Inputs: 
%    dipoles - this may be either the same dipole structure given as input to 
%              the dipplot() function, a 3 by n array of dipole localization or
%              a cell array containing arguments for the dipplot function. Note that
%              the 'coordformat' option below defines the coordinate space for these
%              dipoles (default is MNI). See help dipplot for more information.
%
% Optional 'key', val input pairs:
%  'mri'        - [string or struct] mri file (matlab format or file format read 
%                 by fcdc_read_mri). See dipplot.m help for more information.
%  'method'     - ['alldistance'|'distance'|'entropy'|'relentropy'] method for 
%                            computing density: 
%                 'alldistance' - {default} take into account the gaussian-weighted 
%                            distances from each voxel to all the dipoles. See 
%                            'methodparam' (below) to specify a standard deviation 
%                            (in mm) for the gaussian weight kernel.
%                 'distance' - take into account only the distances to the nearest
%                              dipole for each subject. See 'methodparam' (below).
%                 'entropy' - taking into account only the nearest dipole to each 
%                             voxel for each subject. See 'methodparam' below. 
%                 'relentropy' - as in 'entropy,' but take into account all the 
%                             dipoles for each subject. 
% 'methodparam' - [number] for 'distance'|'alldistance' methods (see above), the
%                 standard deviation (in mm) of the 3-D gaussian smoothing kernel.
%                 For 'entropy'|'relentropy' methods, the number of closest dipoles 
%                 to include {defaults: 20 mm | 20 dipoles }
% 'subsample'   - [integer] subsampling of native MNI image {default: 2 -> 2x2x2}
% 'weight'      - [(1,ncomps) array] for 'distance'|'alldistance' methods, the 
%                 relative weight of each component dipole {default: ones()}
% 'coordformat' - ['mni'|'spherical'] coordinate format if dipole location or 
%                 a structure is given as input. Default is 'mni'.
% 'subjind'     - [(1,ncomps) array] subject index for each dipole model. If two 
%                 dipoles are in one component model, give only one subject index. 
% 'nsessions'   - [integer] for 'alldistance' method, the number of sessions to 
%                 divide the output values by, so that the returned measure is 
%                 dipole density per session {default: 1}
% 'plot'        - ['on'|'off'] force plotting dipole density|entropy 
%                 {default: 'on' if no output arguments, else 'off'}
% 'dipplot'     - ['on'|'off'] plot the dipplot image (used for converting
%                 coordinates (default is 'off')
% 'plotargs'    - {cell array} plotting arguments for mri3dplot() function.
% 'volmesh_fname' - [string] precomputed mesh volume file name. If not
%                 given as input the function will recompute it (it can take from
%                 five to 20 minutes). By default this function save the volume file 
%                 mesh into a file named volmesh_local.mat in the current
%                 folder.
% 'norm2JointProb' - ['on'|'off'] Use joint probability (i.e. sum of all
%                    voxel values == 1) instead of number of dipoles/cm^3.
%                    Should be used for group comparison. (default 'off')
%
% Outputs:
%  dens3d       - [3-D num array] density in dipoles per cubic centimeter. If output
%                 is returned, no plot is produced unless 'plot','on' is specified. 
%  mri          - {MRI structure} used in mri3dplot().
%
% Example: 
%         >> fakedipoles = (rand(3,10)-0.5)*80;
%         >> [dens3d mri] = dipoledensity( fakedipoles, 'coordformat', 'mni'); 
%         >> mri3dplot(dens3d,mri); % replot if no output is given above
%                                   % function is called automatically
%
% ------------------------------------
% NOTES:   to do multiple subject causal-weighted density map, 
% (1) concatenate dipplot coord matrices for all subject
% (2) make g.subjind vector [ones(1,ncompsS1) 2*ones(1,ncompsS2) ... N*ones(1,ncompssN)]
% (3) concatenate normalized outflows for all subjects to form weight vector
% (4) call dipoledensity function with method = 'entropy' or 'relentropy'
% ------------------------------------
%
% See also:
%           EEGLAB: dipplot(), mri3dplot(), Fieldtrip: ft_inside_headmodel() 
%
% Authors: Arnaud Delorme & Scott Makeig SCCN, INC, UCSD
% Modified by: Makoto Miyakoshi
%              Ramon Martinez-Cancino
% Copyright (C) Arnaud Delorme & Scott Makeig, SCCN/INC/UCSD, 2003-
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [prob3d, mri] = dipoledensity(dipplotargs, varargin)

    % TO DO: return in dipplot() the real 3-D location of dipoles (in posxyz)
    %        FIX the dimension order here
    
prob3d = []; mri = [];
if nargin < 1
   help dipoledensity
   return
end

g = finputcheck(varargin, { 'subjind'     'integer'  []               [];
                            'method' 'string' { 'relentropy','entropy','distance','alldistance' } 'alldistance';
                            'methodparam' 'real'     []               20; 
                            'weight'      { 'real','cell' }  []               [];
                            'smooth'      'real'     []               0;
                            'nsessions'   'integer'  []               1;
                            'subsample'   'integer'  []               2;
                            'plotargs'    'cell'     []               {};
                            'plot'        'string'  { 'on','off' }    fastif(nargout == 0, 'on', 'off');
                            'dipplot'     'string'  { 'on','off' }   'off';
                            'coordformat' 'string'  { 'mni','spherical' }   'mni';
                            'normalization' 'string'  { 'on','off' } 'on';
                            'volmesh_fname' 'string'  []  'volmesh_local.mat';
                            'mri'         { 'struct','string' } [] '';
                            'norm2JointProb' 'string'  { 'on','off' } 'off'});
if ischar(g), error(g); end
if ~strcmpi(g.method, 'alldistance') && isempty(g.subjind)
    error('Subject indices are required for this method');
end
if ~iscell(g.weight), g.weight = { g.weight }; end

% plotting dipplot
% ----------------
if ~iscell(dipplotargs) % convert input
    if ~isstruct(dipplotargs)
        if size(dipplotargs,1) == 3, dipplotargs = dipplotargs';
        elseif size(dipplotargs,2) ~= 3
            error('If an array of dipoles is given as entry, there must be 3 columns or 3 rows for x y z');
        end
        model = [];
        for idip = 1:length(dipplotargs)
            model(idip).posxyz = dipplotargs(idip,:);
            model(idip).momxyz = [1 0 0];
            model(idip).rv = 0.5;
        end
        dipplotargs = model;
    end
    dipplotargs = { dipplotargs 'coordformat' g.coordformat };
else 
    dipplotargs = { dipplotargs{:} 'coordformat' g.coordformat };
end
struct = dipplot(dipplotargs{:}, 'plot', g.dipplot, 'density', 'off');
if nargout == 0
    drawnow;
end

% retrieve coordinates in MNI space
% ---------------------------------
if 0 % deprecated
     % find dipoles 
     % ------------
    hmesh = findobj(gcf, 'tag', 'mesh');
    if isempty(hmesh), error('Current figure must contain dipoles'); end
    hh = [];
    disp('Finding dipoles...');
    dips = zeros(1,200);
    for index = 1:1000
        hh = [ hh(:); findobj(gcf, 'tag', ['dipole' int2str(index) ]) ];
        dips(index) = length(findobj(gcf, 'tag', ['dipole' int2str(index) ]));
    end
    
    disp('Retrieving dipole positions ...');
    count = 1;
    for index = 1:length(hh)
        tmp = get(hh(index), 'userdata');
        if length(tmp) == 1
            allx(count) = tmp.eleccoord(1,1);
            ally(count) = tmp.eleccoord(1,2);
            allz(count) = tmp.eleccoord(1,3);
            alli(count) = index;
            count = count + 1;
        end
    end
end    

% check weights
% -------------
if ~isempty(g.weight{1})
    if ~iscell(g.weight)
        if length(g.weight) ~= length(struct)
            error('There must be as many elements in the weight matrix as there are dipoles')
        end
    else
        if length(g.weight{1}) ~= length(struct) || length(g.weight{1}) ~= length(g.weight{end})
            error('There must be as many elements in the weight matrix as there are dipoles')
        end
    end
else
    g.weight = { ones( 1, length(struct)) };
end
if ~isempty(g.subjind)
    if length(g.subjind) ~= length(struct)
        error('There must be as many element in the subject matrix as there are dipoles')
    end
else
    g.subjind = ones( 1, length(struct));
end

% decoding dipole locations
% -------------------------
disp('Retrieving dipole positions ...');
count = 1;
for index = 1:length(struct)
    dips = size(struct(index).eleccoord,1);
    for dip = 1:dips
        allx(count) = struct(index).eleccoord(dip,1);
        ally(count) = struct(index).eleccoord(dip,2);
        allz(count) = struct(index).eleccoord(dip,3);
        alli(count) = index;
        allw1(count) = g.weight{1}(  index)/dips;
        allw2(count) = g.weight{end}(index)/dips;
        alls(count) = g.subjind(index);
        count = count + 1;
    end
end
g.weight{1}    = allw1;
g.weight{end}  = allw2;
g.subjind = alls;

% read MRI file
% -------------
if isempty(g.mri) % default MRI file
    dipfitdefs;
    load('-mat', template_models(1).mrifile); % load mri variable
    g.mri = mri;
end
if ischar(g.mri)
    try
        mri = load('-mat', g.mri);
        mri = mri.mri;
    catch
        disp('Failed to read Matlab file. Attempt to read MRI file using function read_fcdc_mri');
        try
            warning off;
            mri = read_fcdc_mri(g.mri);
            mri.anatomy = round(gammacorrection( mri.anatomy, 0.8));
            mri.anatomy = uint8(round(mri.anatomy/max(reshape(mri.anatomy, prod(mri.dim),1))*255));
            % WARNING: if using double instead of int8, the scaling is different 
            % [-128 to 128 and 0 is not good]
            % WARNING: the transform matrix is not 1, 1, 1 on the diagonal, some slices may be 
            % misplaced
            warning on;
        catch
            error('Cannot load file using read_fcdc_mri');
        end
    end
    g.mri = mri; % output the anatomic mri image 
end


% reserve array for density
% -------------------------
prob3d = {zeros(ceil(g.mri.dim/g.subsample)) };
for i = 2:length(g.weight), prob3d{i} = prob3d{1}; end

% compute voxel size
% ------------------
point1 = g.mri.transform * [ 1 1 1 1 ]';
point2 = g.mri.transform * [ 2 2 2 1 ]';
voxvol = sum((point1(1:3)-point2(1:3)).^2)*g.subsample^3; % in mm

% compute global subject entropy if necessary
% -------------------------------------------
vals   = unique_bc(g.subjind); % the unique subject indices
if strcmpi(g.method, 'relentropy') || strcmpi(g.method, 'entropy') %%%%% entropy %%%%%%%
    newind = zeros(size(g.subjind));
    for index = 1:length(vals) % foreach subject in the cluster
        tmpind = find(g.subjind == vals(index)); % dipoles for the subject
        totcount(index) = length(tmpind); % store the number of subject dipoles
        newind(tmpind) = index; % put subject index into newind
    end
    g.subjind = newind;
    gp = totcount/sum(totcount);
    globent = -sum(gp.*log(gp));
end

% compute volume inside head mesh
% -------------------------------
dipfitdefs; % get the location of standard BEM volume file
tmp = load('-mat',DIPOLEDENSITY_STDBEM); % load MNI mesh

if ~exist(g.volmesh_fname)
    % Checking for Fieldtrip
    if exist('ft_electroderealign', 'file')~=2,error('dipoledensity: Fieldtrip toolbox is required'); end
    
    disp('Computing volume within head mesh...');
    [X Y Z]           = meshgrid(g.mri.xgrid(1:g.subsample:end)+g.subsample/2, ...
                                 g.mri.ygrid(1:g.subsample:end)+g.subsample/2, ...
                                 g.mri.zgrid(1:g.subsample:end)+g.subsample/2);
    [indX indY indZ ] = meshgrid(1:length(g.mri.xgrid(1:g.subsample:end)), ...
                                 1:length(g.mri.ygrid(1:g.subsample:end)), ...
                                 1:length(g.mri.zgrid(1:g.subsample:end)));
    allpoints = [ X(:)'    ; Y(:)'   ; Z(:)' ];
    allinds   = [ indX(:)' ; indY(:)'; indZ(:)' ];
    allpoints = g.mri.transform * [ allpoints ; ones(1, size(allpoints,2)) ];
    allpoints(4,:) = [];

    olddir = pwd;
    tmppath = which('ft_electroderealign');
    tmppath = fullfile(fileparts(tmppath), 'private');
    cd(tmppath);
    inside = ft_inside_headmodel(allpoints', tmp.vol);
    Inside = find(inside); Outside = find(~inside);
    cd(olddir);
    disp('Done.');
    
    if 0 % old code using Delaunay %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        P = tmp.vol.bnd(1).pnt;
        T = delaunayn(P); % recompute triangularization (the original one is not compatible 
                          % with tsearchn) get coordinates of all points in the volume
        % search for points inside or outside the volume (takes about 14 minutes!)
        IO = tsearchn(P, T, allpoints');
        Inside        = find(isnan(IO));
        Outside       = find(~isnan(IO));
        disp('Done.');
    end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try 
        save('-mat', g.volmesh_fname, 'allpoints', 'allinds', 'Inside', 'Outside');
        disp('Saving file containing inside/outide voxel indices...');
    catch, end
else
    disp('Loading file containing inside/outide voxel indices...');
    load('-mat',g.volmesh_fname);
end
InsidePoints  = allpoints(:, Inside);
InsideIndices = allinds(:, Inside);

% scan grid and compute entropy at each voxel
% -------------------------------------------
edges = [0.5:1:length(vals)+0.5];

if ~strcmpi(g.method, 'alldistance') 
    fprintf('Computing (of %d):', size(InsideIndices,2));
    % entropy calculation: have to scan voxels
    % ----------------------------------------
    for i = 1:size(InsideIndices,2)
        
        alldists = (InsidePoints(1,i) - allx).^2 ...
                 + (InsidePoints(2,i) - ally).^2 ...
                 + (InsidePoints(3,i) - allz).^2;
        [tmpsort indsort] = sort(alldists); % sort dipoles by distance
        tmpweights{1}   = g.weight{1}(  indsort);
        tmpweights{end} = g.weight{end}(indsort);
       
        if strcmpi(g.method, 'relentropy') || strcmpi(g.method, 'entropy') %%%%% entropy %%%%%%%
            
            subjs  = g.subjind(indsort(1:g.methodparam)); % get subject indices of closest dipoles
            p      = histc(subjs, edges);
            if strcmpi(g.method, 'relentropy')
                p      = p(1:end-1)./totcount; 
                % this should be uniform if p conforms to global count for all subjects
            end
            p      = p/sum(p);
            p(find(p == 0)) = [];
            for tmpi = 1:length(g.weight)
                prob3d{1}(InsideIndices(1,i), InsideIndices(2,i), InsideIndices(3,i)) = -sum(p.*log(p));
            end
        else
            % distance to each subject
            ordsubjs  = g.subjind(indsort);
            for index = 1:length(vals) % for each subject
                tmpind = find(ordsubjs == vals(index));
                if strcmpi(g.method,'distance')
                    use_dipoles(index) = tmpind(1); % find their nearest dipole 
                end
            end
            for tmpi = 1:length(g.weight)
                prob3d{tmpi}(InsideIndices(1,i), InsideIndices(2,i), InsideIndices(3,i)) = ...
                   sum(tmpweights{tmpi}(use_dipoles).*exp(-tmpsort(use_dipoles)/ ...
                           (2*g.methodparam^2))); % 3-D gaussian smooth
            end
        end
        if mod(i,100) == 0, fprintf('%d ', i); end
    end
else % 'alldistance'
    % distance calculation: can scan dipoles instead of voxels (since linear)
    % --------------------------------------------------------
    %alldists = allx.^2 + ally.^2 + allz.^2;
    %figure; hist(alldists); return; % look at distribution of distances
    
    fprintf('Computing (of %d):', size(allx,2));
    for tmpi=1:length(g.weight)
        tmpprob{tmpi} = zeros(1, size(InsidePoints,2));
    end
    if length(g.weight) > 1, tmpprob2 = tmpprob; end
    for i = 1:size(allx,2)
        alldists = (InsidePoints(1,:) - allx(i)).^2 + ...
                   (InsidePoints(2,:) - ally(i)).^2 + ...
                   (InsidePoints(3,:) - allz(i)).^2;
%         alldists = 1;  % TM
        for tmpi=1:length(g.weight)
            tmpprob{tmpi} = tmpprob{tmpi} + g.weight{tmpi}(i)*exp(-alldists/(2*g.methodparam^2)); % 3-D gaussian smooth
            if any(isinf(tmpprob{tmpi})), error('Infinite value in probability calculation'); end
        end
        if mod(i,50) == 0, fprintf('%d ', i); end
    end
    % copy values to 3-D mesh
    % -----------------------
    for i = 1:length(Inside)
        pnts = allinds(:,Inside(i));
        for tmpi = 1:length(g.weight)
            prob3d{tmpi}(pnts(1), pnts(2), pnts(3)) = tmpprob{tmpi}(i);
        end
    end
    
end
fprintf('\n');

% normalize for points inside and outside the volume
% norm2JointProb is applied before plotting
% --------------------------------------------------
if strcmpi(g.method, 'alldistance') && strcmpi(g.normalization,'on')
    for i =1:length(g.weight)
        disp('Normalizing to dipole/mm^3');
        if any(prob3d{i}(:)<0)
            fprintf('WARNING: Some probabilities are negative, this will likely cause problems when normalizing probabilities.\n');
            fprintf('It is highly recommended to turn normaliziation off by using ''normalization'' key to ''off''.\n');
        end
        totval = sum(prob3d{i}(:));                                     % total values in the head
        totdip = size(allx,2);                                          % number of dipoles
        prob3d{i} = (prob3d{i}/totval*totdip/voxvol*1000)/g.nsessions;  % time 1000 to get cubic centimeters
    end
end

% resample matrix
% ----------------
if g.subsample ~= 1
    for i =1:length(g.weight)
        prob3d{i} = prob3d{i}/g.subsample;
        newprob3d = zeros(g.mri.dim);
        X = ceil(g.mri.xgrid/g.subsample);
        Y = ceil(g.mri.ygrid/g.subsample);
        Z = ceil(g.mri.zgrid/g.subsample);
        for index = 1:size(newprob3d,3)
            newprob3d(:,:,index) = prob3d{i}(X,Y,Z(index));
        end    
        prob3d{i} = newprob3d;
    end
end

% 3-D smoothing
% -------------
if g.smooth ~= 0
    disp('Smoothing...');
    for i =1:length(g.weight)
        prob3d{i} = smooth3d(prob3d{i}, g.smooth);
    end
end

% Perform normalization so that the total sum of joint prob == 1
if strcmpi(g.norm2JointProb, 'on')
    prob3d{i} = prob3d{i}/sum(prob3d{i}(:));
end

% plotting
% --------
if strcmpi(g.plot, 'off')
    close gcf;
else
    mri3dplot( prob3d, g.mri, g.plotargs{:}); % plot the density using mri3dplot()
end
return;
