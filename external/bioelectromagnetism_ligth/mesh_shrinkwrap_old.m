0001 function [FV, Edges] = mesh_shrinkwrap(vol,FV,smooth,vthresh,interpVal,...
0002     fitval,fittol,fititer,fitchange,fitvattr)
0003 
0004 % mesh_shrinkwrap - Tesselate the surface of a 3D volume
0005 %
0006 % [FV, Edges] = mesh_shrinkwrap(vol,FV,smooth,vthresh,interpVal,...
0007 %               fitval,fittol,fititer,fitchange,fitvattr)
0008 %
0009 % vol       - a 3D image volume
0010 % FV        - input tesselation; if empty, sphere tesselation
0011 %             is created.  FV has fields FV.vertices, FV.faces
0012 % smooth    - Gaussian smoothing of vol (5x5x5 kernel), 0|1 (default 0)
0013 % vthresh   - binarise vol at fitval, 0|1 (default 0)
0014 % interpVal - use radial interpolation (faster) or incremental
0015 %             radial shrink (slower), 0|1 (default 1, faster)
0016 % fitval    - image intensity to shrink wrap (default 20)
0017 % fittol    - image intensity tolerance (default 5)
0018 % fititer   - max number of iterations to fit (default 200)
0019 % fitchange - least sig. change in intensity values
0020 %             between fit iterations (default 2)
0021 % fitvattr  - vertex attraction (constraint), 0:1, smaller
0022 %             values are less constraint; close to 0 for
0023 %             no constraint is useful when dealing with
0024 %             binary volumes, otherwise 0.4 (40%) seems OK
0025 %
0026 % FV        - a struct with 2562 vertices and 5120 faces
0027 % Edges     - a [2562,2562] matrix of edge connectivity for FV
0028 %
0029 % An alternative to isosurface for volumes with an external
0030 % surface that can be shrink-wrapped. It has been developed to
0031 % find the scalp surface for MRI of the human head.
0032 %
0033 % It starts with a sphere tesselation (large radius) and moves
0034 % each vertex point toward the centre of the volume until it
0035 % lies at or near the fitval.  The function is not optimised
0036 % for speed, but it should produce reasonable results.
0037 %
0038 % Example of creating a scalp tesselation for SPM T1 MRI template:
0039 %
0040 %   avw = avw_read('T1');
0041 %   FV = mesh_shrinkwrap(avw.img,[],0,0,intensity,5.0,50,0.5,0.4);
0042 %   patch('vertices',FV.vertices,'faces',FV.faces,'facecolor',[.6 .6 .6]);
0043 %
0044 % Example of creating a skull from FSL BET skull volume:
0045 %
0046 %   avw = avw_read('T1_skull');
0047 %   FV = mesh_shrinkwrap(avw.img,[],1,0,intensity,0.2,10,0.005,0.1);
0048 %   patch('vertices',FV.vertices,'faces',FV.faces,'facecolor',[.6 .6 .6]);
0049 %
0050 % ***** NOTE *****
0051 % An important limitation at present is that it doesn't
0052 % read the image dimensions, but assumes they are 1mm^3.  Further
0053 % versions might take a standard analyze volume and check the
0054 % header to scale the result according to the image dimensions. At
0055 % present, this must be done with the results of this function.  The
0056 % output vertex coordinates are in mm with an origin at (0,0,0),
0057 % which lies at the center of the MRI volume.
0058 %
0059 % See also: ISOSURFACE, SPHERE_TRI, MESH_REFINE_TRI4,
0060 %           MESH_BEM_SHELLS_FUNC, MESH_BEM_SHELLS_SCRIPT
0061 %
0062 
0063 % $Revision: 1.1 $ $Date: 2009-04-28 22:13:58 $
0064 
0065 % Licence:  GNU GPL, no implied or express warranties
0066 % History:  07/2002, Darren.Weber_at_radiology.ucsf.edu
0067 %                    - created to provide alternative to
0068 %                      isosurface in matlab R12
0069 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
0070 
0071 version = '[$Revision: 1.1 $]';
0072 fprintf('MESH_SHRINKWRAP [v%s]\n',version(12:16));  tic;
0073 
0074 % Parse arguments
0075 
0076 if ~exist('vol','var'), error('...no input volume\n');
0077 elseif isempty(vol),    error('...empty input volume\n');
0078 end
0079 if isstruct(vol),
0080     if isfield(vol,'img'), vol = vol.img;
0081     else
0082         error('...input volume is a struct, it must be a 3D image volume only\n');
0083     end
0084 end
0085 
0086 if ~exist('smooth','var'),   smooth = 0;
0087 elseif isempty(smooth),      smooth = 0;
0088 end
0089 
0090 if ~exist('vthresh','var'),  vthresh = 0;
0091 elseif isempty(vthresh),     vthresh = 0;
0092 end
0093 
0094 if ~exist('interpVal','var'), interpVal = 1;
0095 elseif isempty(interpVal),    interpVal = 1;
0096 end
0097 
0098 if ~exist('fitval','var'),   fit.val = 20;
0099 elseif isempty(fitval),      fit.val = 20;
0100 else                         fit.val = fitval;
0101 end
0102 
0103 if ~exist('fittol','var'),   fit.tol = 5;
0104 elseif isempty(fittol),      fit.tol = 5;
0105 else                         fit.tol = fittol;
0106 end
0107 
0108 if fit.val <= fit.tol,
0109     error('...must use fit tolerance < fit value\n');
0110 end
0111 
0112 if ~exist('fititer','var'),  fit.iter = 200;
0113 elseif isempty(fititer),     fit.iter = 200;
0114 else                         fit.iter = fititer;
0115 end
0116 
0117 if ~exist('fitchange','var'),fit.change = 2;
0118 elseif isempty(fitchange),   fit.change = 2;
0119 else                         fit.change = fitchange;
0120 end
0121 
0122 if ~exist('fitvattr','var'), fit.vattr = 0.4;
0123 elseif isempty(fitvattr),    fit.vattr = 0.4;
0124 else                         fit.vattr = fitvattr;
0125 end
0126 if fit.vattr > 1,
0127     fprintf('...fit vertattr (v) must be 0 <= v <= 1, setting v = 1\n');
0128     fit.vattr = 1;
0129 end
0130 if fit.vattr < 0,
0131     fprintf('...fit vertattr (v) must be 0 <= v <= 1, setting v = 0.\n');
0132     fit.vattr = 0;
0133 end
0134 
0135 % MAIN
0136 
0137 % Find approximate centre of volume
0138 xdim = size(vol,1);
0139 ydim = size(vol,2);
0140 zdim = size(vol,3);
0141 
0142 origin(1) = floor(xdim/2);
0143 origin(2) = floor(ydim/2);
0144 origin(3) = floor(zdim/2);
0145 
0146 % Check whether to create a sphere tesselation
0147 % or use an input tesselation as the start point
0148 sphere = 0;
0149 if ~exist('FV','var'),
0150     sphere = 1;
0151 elseif ~isfield(FV,'vertices'),
0152     sphere = 1;
0153 elseif ~isfield(FV,'faces'),
0154     sphere = 1;
0155 elseif isempty(FV.vertices),
0156     sphere = 1;
0157 elseif isempty(FV.faces),
0158     sphere = 1;
0159 end
0160 if sphere,
0161     % Create a sphere tesselation to encompass the volume
0162     radius = max([xdim ydim zdim]) / 1.5;
0163     FV = sphere_tri('ico',4,radius); % 2562 vertices
0164 else
0165     fprintf('...using input FV tesselation...\n');
0166 end
0167 
0168 
0169 % the 'edge' matrix is the connectivity of all vertices,
0170 % used to find neighbours during movement of vertices,
0171 % including smoothing the tesselation
0172 FV.edge = mesh_edges(FV);
0173 
0174 
0175 % Shift the centre of the sphere to the centre of the volume
0176 centre = repmat(origin,size(FV.vertices,1),1);
0177 FV.vertices = FV.vertices + centre;
0178 
0179 
0180 % 'Binarise' the volume, removing all values below
0181 % a threshold, setting all others to threshold
0182 if vthresh,
0183     fprintf('...thresholding volume...'); tic;
0184     Vindex = find(vol < fit.val);
0185     vol(Vindex) = 0;
0186     Vindex = find(vol >= fit.val);
0187     vol(Vindex) = fit.val;
0188     t = toc; fprintf('done (%5.2f sec)\n',t);
0189 end
0190 binvol = find(vol > 1);
0191 
0192 % smooth the volume
0193 if smooth,
0194     fprintf('...gaussian smoothing (5-10 minutes)...'); tic;
0195     vol = smooth3(vol,'gaussian',5,.8);
0196     t = toc; fprintf('done (%5.2f sec)\n',t);
0197 end
0198 
0199 
0200 % Now begin recursion
0201 fprintf('...fitting...\n');    tic;
0202 
0203 i = 1;
0204 Fminima = 0;
0205 intensityAtDMean = [0 0];
0206 
0207 while i <= fit.iter,
0208     
0209     if interpVal,
0210         % use radial interpolation method, moving directly
0211         % to the intensity value nearest correct intensity
0212         [FV, intensityAtD, D] = locate_val(FV,vol,origin,fit);
0213     else
0214         % use incremental method, moving along radial line
0215         % gradually until finding correct intensity
0216         [FV, intensityAtD, D] = shrink_wrap(FV,vol,origin,fit);
0217     end
0218     
0219     intensityAtDMean(1) = intensityAtDMean(2);
0220     intensityAtDMean(2) = mean(intensityAtD);
0221     
0222     fprintf('...distance:  mean = %8.4f mm, std = %8.4f mm\n',mean(D),std(D));
0223     fprintf('...intensity: mean = %8.4f,    std = %8.4f\n',...
0224         mean(intensityAtD),std(intensityAtD));
0225     fprintf('...real iteration: %3d\n',i);
0226     
0227     % Is the mean distance reasonable?
0228     if mean(D) < 0.5,
0229         error('...mean distance < 0.5 mm!\n');
0230     end
0231     
0232     % MDifVal is the mean of the absolute difference
0233     % between the vertex intensity and the fit intensity
0234     MDifVal = abs(intensityAtDMean(2) - fit.val);
0235     
0236     % Is the mean difference within the tolerance range?
0237     if MDifVal < fit.tol,
0238         fprintf('...mean intensity difference < tolerance (%5.2f +/- %5.2f)\n',...
0239             fit.val,fit.tol);
0240         break;
0241     else
0242         fprintf('...mean intensity difference > tolerance (%5.2f +/- %5.2f)\n',...
0243             fit.val,fit.tol);
0244     end
0245     
0246     % How much has the intensity values changed?
0247     if (i > 1) & intensityAtDMean(2) > 0,
0248         if intensityAtDMean(2) - intensityAtDMean(1) < fit.change,
0249             fprintf('...no significant intensity change (< %5.2f) in this iteration\n',...
0250                 fit.change);
0251             Fminima = Fminima + 1;
0252             if Fminima >= 5,
0253                 fprintf('...no significant intensity change in last 5 iterations\n');
0254                 break;
0255             end
0256         else
0257             Fminima = 0;
0258         end
0259     end
0260     
0261     % Ensure that iterations begin when MDifVal is truly initialised
0262     if isnan(MDifVal),
0263         i = 1;
0264     else,
0265         i = i + 1;
0266     end
0267     
0268 end
0269 
0270 FV = mesh_smooth(FV,origin,fit.vattr);
0271 
0272 % Remove large edges matrix from FV
0273 Edges = FV.edge;
0274 FV = struct('vertices',FV.vertices,'faces',FV.faces);
0275 
0276 % Now center the output vertices at 0,0,0 by subtracting
0277 % the volume centroid
0278 FV.vertices(:,1) = FV.vertices(:,1) - origin(1);
0279 FV.vertices(:,2) = FV.vertices(:,2) - origin(2);
0280 FV.vertices(:,3) = FV.vertices(:,3) - origin(3);
0281 
0282 t=toc; fprintf('...done (%5.2f sec).\n\n',t);
0283 
0284 return
0285 
0286 
0287 
0288 
0289 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
0290 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
0291 function [FV, intensityAtD, D] = locate_val(FV,vol,origin,fit),
0292     
0293     xo = origin(1); yo = origin(2); zo = origin(3);
0294     
0295     Nvert = size(FV.vertices,1);
0296     progress = round(.1 * Nvert);
0297     
0298     % Initialise difference intensity & distance arrays
0299     intensityAtD = zeros(Nvert,1);
0300     D    = intensityAtD;
0301     
0302     % Find distance and direction cosines for line from
0303     % origin to all vertices
0304     XV = FV.vertices(:,1);
0305     YV = FV.vertices(:,2);
0306     ZV = FV.vertices(:,3);
0307     DV = sqrt( (XV-xo).^2 + (YV-yo).^2 + (ZV-zo).^2 );
0308     LV = (XV-xo)./DV; % cos alpha
0309     MV = (YV-yo)./DV; % cos beta
0310     NV = (ZV-zo)./DV; % cos gamma
0311     
0312     % Check for binary volume data, if empty, binary
0313     binvol = find(vol > 1);
0314     
0315     % Locate each vertex at a given fit value
0316     tic
0317     for v = 1:Nvert,
0318         
0319         if v > progress,
0320             fprintf('...interp3 processed %4d of %4d vertices',progress,Nvert);
0321             t = toc; fprintf(' (%5.2f sec)\n',t);
0322             progress = progress + progress;
0323         end
0324         
0325         % Find direction cosines for line from origin to vertex
0326         x = XV(v);
0327         y = YV(v);
0328         z = ZV(v);
0329         d = DV(v);
0330         l = LV(v); % cos alpha
0331         m = MV(v); % cos beta
0332         n = NV(v); % cos gamma
0333         
0334         % find discrete points between origin
0335         % and vertex + 50% of vertex distance
0336         points = 250;
0337         
0338         Darray = linspace(0,(d + .2 * d),points);
0339         
0340         L = repmat(l,1,points);
0341         M = repmat(m,1,points);
0342         N = repmat(n,1,points);
0343         
0344         XI = (L .* Darray) + xo;
0345         YI = (M .* Darray) + yo;
0346         ZI = (N .* Darray) + zo;
0347         
0348         % interpolate volume values at these points
0349         VI = interp3(vol,YI,XI,ZI,'*linear');
0350         
0351         % do we have a binary volume (no values > 1)
0352         if isempty(binvol),
0353             maxindex = max(find(VI>0));
0354             if maxindex,
0355                 D(v) = Darray(maxindex);
0356             end
0357         else
0358             % find the finite values of VI
0359             index = max(find(VI(isfinite(VI))));
0360             if index,
0361                 
0362                 % Find nearest volume value to the required fit value
0363                 nearest = max(find(VI >= fit.val));
0364                 
0365                 %[ nearest, value ] = NearestArrayPoint( VI, fit.val );
0366                 
0367                 % Check this nearest index against a differential
0368                 % negative peak value
0369                 %diffVI = diff(VI);
0370                 %if max(VI) > 1,
0371                 %    diffindex = find(diffVI < -20);
0372                 %else
0373                 % probably a binary volume
0374                 %    diffindex = find(diffVI < 0);
0375                 %end
0376                 
0377                 % now set d
0378                 if nearest,
0379                     D(v) = Darray(nearest);
0380                 end
0381             end
0382         end
0383         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
0384         % Constrain relocation by fit.vattr,
0385         % some % of distance from neighbours
0386         
0387         vi = find(FV.edge(v,:));  % the neighbours' indices
0388         X = FV.vertices(vi,1);    % the neighbours' vertices
0389         Y = FV.vertices(vi,2);
0390         Z = FV.vertices(vi,3);
0391         
0392         % Find neighbour distances
0393         DN = sqrt( (X-xo).^2 + (Y-yo).^2 + (Z-zo).^2 );
0394         % Find mean distance of neighbours
0395         DNmean = mean(DN);
0396         
0397         minattr = fit.vattr;
0398         maxattr = 1 + fit.vattr;
0399         
0400         if D(v) < (minattr * DNmean),
0401             D(v) = minattr * DNmean;
0402         end
0403         if D(v) > (maxattr * DNmean),
0404             D(v) = maxattr * DNmean;
0405         end
0406         if D(v) == 0, D(v) = DNmean; end
0407         
0408         % relocate vertex to new distance
0409         x = (l * D(v)) + xo;
0410         y = (m * D(v)) + yo;
0411         z = (n * D(v)) + zo;
0412         
0413         FV.vertices(v,:) = [ x y z ];
0414         
0415         % Find intensity value at this distance
0416         intensityAtD(v) = interp1(Darray,VI,D(v),'linear');
0417         
0418     end
0419     
0420     if isempty(binvol),
0421         % check outliers and smooth twice for binary volumes
0422         FV = vertex_outliers(FV, D, origin);
0423         FV = mesh_smooth(FV,origin,fit.vattr);
0424     end
0425     FV = mesh_smooth(FV,origin,fit.vattr);
0426     
0427 return
0428 
0429 
0430 
0431 
0432 
0433 
0434 
0435 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
0436 function [FV, intensityAtD, D] = shrink_wrap(FV,vol,origin,fit),
0437     
0438     xo = origin(1); yo = origin(2); zo = origin(3);
0439     
0440     Nvert = size(FV.vertices,1);
0441     
0442     intensityAtD = zeros(Nvert,1);
0443     D    = intensityAtD;
0444     
0445     for v = 1:Nvert,
0446         
0447         x = FV.vertices(v,1);
0448         y = FV.vertices(v,2);
0449         z = FV.vertices(v,3);
0450         
0451         % Find distance of vertex from origin
0452         d = sqrt( (x-xo)^2 + (y-yo)^2 + (z-zo)^2 );
0453         
0454         % check whether vertex is already at fit.val in vol
0455         
0456         volval = vol_val(vol,x,y,z);
0457         
0458         if isnan(volval) | (volval < fit.val - fit.tol) | (volval > fit.val + fit.tol),
0459             
0460             % Find direction cosines for line from centre to vertex
0461             l = (x-xo)/d; % cos alpha
0462             m = (y-yo)/d; % cos beta
0463             n = (z-zo)/d; % cos gamma
0464             
0465             % now modify d by fit.dist
0466             if isnan(volval) | (volval < fit.val - fit.tol),
0467                 d = d - (d * fit.vattr);
0468             else
0469                 d = d + (d * fit.vattr);
0470             end
0471             
0472             
0473             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
0474             % Constrain relocation by fit.vattr,
0475             % some % of distance from neighbours
0476             
0477             vi = find(FV.edge(v,:));  % the neighbours' indices
0478             X = FV.vertices(vi,1);    % the neighbours' vertices
0479             Y = FV.vertices(vi,2);
0480             Z = FV.vertices(vi,3);
0481             
0482             % Find neighbour distances
0483             DN = sqrt( (X-xo).^2 + (Y-yo).^2 + (Z-zo).^2 );
0484             % Find mean distance of neighbours
0485             DNmean = mean(DN);
0486             
0487             minattr = fit.vattr;
0488             maxattr = 1 + fit.vattr;
0489             
0490             if d < (minattr * DNmean),
0491                 d = minattr * DNmean;
0492             elseif d > (maxattr * DNmean),
0493                 d = maxattr * DNmean;
0494             end
0495             
0496             % locate vertex at this new distance
0497             x = (l * d) + xo;
0498             y = (m * d) + yo;
0499             z = (n * d) + zo;
0500             
0501             FV.vertices(v,:) = [ x y z ];
0502         end
0503         
0504         intensityAtD(v) = vol_val(vol,x,y,z);
0505         D(v) = d;
0506     end
0507     
0508     FV = mesh_smooth(FV,origin,fit.vattr);
0509     
0510 return
0511 
0512 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
0513 % Currently not calling this function (Oct 02)
0514 function [val] = vol_val(vol,x,y,z),
0515     
0516     % This function just ensures that xyz are
0517     % actually within the volume before trying
0518     % to get a volume value
0519     
0520     val = nan; % assume zero value
0521     
0522     x = round(x);
0523     y = round(y);
0524     z = round(z);
0525     
0526     if x > 0 & y > 0 & z > 0,
0527         
0528         % get bounds of vol
0529         Xv = size(vol,1);
0530         Yv = size(vol,2);
0531         Zv = size(vol,3);
0532         
0533         if x <= Xv & y <= Yv & z <= Zv,
0534             % OK return volume value at xyz
0535             val = vol(x,y,z);
0536         end
0537     end
0538     
0539 return
0540 
0541 
0542 
0543 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
0544 function [FV] = vertex_outliers(FV, D, origin),
0545     
0546     xo = origin(1); yo = origin(2); zo = origin(3);
0547     
0548     % Screen FV for outlying vertices, using
0549     % mean +/- 2 * stdev of distance from origin
0550     DistMean = mean(D);
0551     DistStDev = std(D);
0552     DistMax = DistMean + 2 * DistStDev;
0553     DistMin = DistMean - 2 * DistStDev;
0554     
0555     for v = 1:size(FV.vertices,1),
0556         
0557         if D(v) >= DistMax,
0558             D(v) = DistMean;
0559             relocate = 1;
0560         elseif D(v) <= DistMin,
0561             D(v) = DistMean;
0562             relocate = 1;
0563         else
0564             relocate = 0;
0565         end
0566         
0567         if relocate,
0568             x = FV.vertices(v,1);
0569             y = FV.vertices(v,2);
0570             z = FV.vertices(v,3);
0571             
0572             % Find direction cosines for line from centre to vertex
0573             d = sqrt( (x-xo)^2 + (y-yo)^2 + (z-zo)^2 );
0574             l = (x-xo)/d; % cos alpha
0575             m = (y-yo)/d; % cos beta
0576             n = (z-zo)/d; % cos gamma
0577             
0578             % relocate vertex to this new distance
0579             x = (l * D(v)) + xo;
0580             y = (m * D(v)) + yo;
0581             z = (n * D(v)) + zo;
0582             
0583             FV.vertices(v,:) = [ x y z ];
0584         end
0585     end
0586 return
