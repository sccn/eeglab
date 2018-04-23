function M2T = mni2tal_matrix()

% MNI2TAL_MATRIX - Talairach to MNI coordinates (best guess)
% 
% MNI2TALAIRACH = mni2tal_matrix
% 
% MNI2TALAIRACH is a struct containing rotation matrices 
% used by mni2tal and tal2mni
% 
% See also, MNI2TAL, TAL2MNI &
% http://www.mrc-cbu.cam.ac.uk/Imaging/mnispace.html
% 

% $Revision: 1.1 $ $Date: 2009-01-30 02:48:33 $

% Licence:  GNU GPL, no express or implied warranties
% Matthew Brett 2/2/01, matthew.brett@mrc-cbu.cam.ac.uk
% modified 02/2003, Darren.Weber_at_radiology.ucsf.edu
%                   - removed dependence on spm_matrix by
%                     creating this function, thereby 
%                     abstracting the important matrix
%                     transforms (easier to change if needed).
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% See notes below for explanations...


% rotn  = spm_matrix([0 0 0 0.05]); % similar to Rx(eye(3),-0.05), DLW
M2T.rotn  = [      1         0         0         0;
                   0    0.9988    0.0500         0;
                   0   -0.0500    0.9988         0;
                   0         0         0    1.0000 ];


% upz   = spm_matrix([0 0 0 0 0 0 0.99 0.97 0.92]);
M2T.upZ   = [ 0.9900         0         0         0;
                   0    0.9700         0         0;
                   0         0    0.9200         0;
                   0         0         0    1.0000 ];


% downz = spm_matrix([0 0 0 0 0 0 0.99 0.97 0.84]);
M2T.downZ = [ 0.9900         0         0         0;
                   0    0.9700         0         0;
                   0         0    0.8400         0;
                   0         0         0    1.0000 ];

% from original mni2tal...
%upT   = spm_matrix([0 0 0 0.05 0 0 0.99 0.97 0.92]);
%downT = spm_matrix([0 0 0 0.05 0 0 0.99 0.97 0.84]);

               

return


% from http://www.mrc-cbu.cam.ac.uk/Imaging/mnispace.html
% 
% Approach 2: a non-linear transform of MNI to Talairach
% 
% An alternative is to use some sort of transformation that
% may differ for different brain areas. One method might be 
% to do an automated non-linear match of the MNI to the 
% Talairach brain. For example, you could apply an SPM or 
% AIR warping algorithm. However, there are two problems 
% here. First, as we stated above, we do not have an MRI 
% image of the brain in the Talairach atlas, which was a 
% post-mortem specimen. Second, the automated non-linear 
% transforms produce quite complex equations relating the 
% two sets of coordinates. 
% 
% An alternative is to apply something like the transform 
% that Talairach and Tournoux designed; here different 
% linear transforms are applied to different brain regions. 
% This is the approach I describe below. 
% 
% To get a good match for both the temporal lobes and the 
% top of the brain, I used different zooms, in the Z (down/up) 
% direction, for the brain above the level of the AC/PC line, 
% and the brain below. The algorithm was:
% 
% I assumed that the AC was in the correct position in the MNI 
% brain, and therefore that no translations were necessary; 
% Assumed that the MNI brain was in the correct orientation in 
% terms of rotation around the Y axis (roll) and the Z axis (yaw); 
% Using the SPM99b display tool, I compared the MNI brain to the 
% images in the Talairach atlas; 
% 
% Compared to the atlas, the MNI brain seemed tipped backwards, 
% so that the cerebellar / cerebral cortex line in the sagittal 
% view, at the AC, was too low. Similarly, the bottom of the 
% anterior part of the corpus collosum seemed too high. I 
% therefore applied a small (0.05 radian) pitch correction to 
% the MNI brain; 
% 
% Matching the top of the MNI brain to the top of the brain in 
% the atlas, required a zoom of 0.92 in Z. Similarly a Y zoom 
% of 0.97 was required as a best compromise in matching the front 
% and back of the MNI brain to the atlas. The left / right match 
% required a 0.99 zoom in X;
% 
% The transform above provided a good match for the brain superior 
% to the AC/PC line, but a poor match below, with the temporal lobes 
% extending further downwards in the MNI brain than in the atlas. I 
% therefore derived a transform for the brain below the AC/PC line, 
% that was the same as the transform above, except with a Z zoom of 
% 0.84; 
% 
% This algorithm gave me the following transformations: 
% 
% Above the AC (Z >= 0): 
% 
% X'= 0.9900X 
% 
% Y'=  0.9688Y +0.0460Z 
% 
% Z'= -0.0485Y +0.9189Z 
% 
% 
% Below the AC (Z < 0): 
% 
% X'= 0.9900X 
% 
% Y'=  0.9688Y +0.0420Z 
% 
% Z'= -0.0485Y +0.8390Z 
% 
% 
% The matlab function mni2tal.m implements these transforms. 
% It returns estimated Talairach coordinates, from the 
% transformations above, for given points in the MNI brain. 
% To use it, save as mni2tal.m somewhere on your matlab path. 
% 
% So, taking our example point in the MNI brain, X = 10mm, Y = 12mm, Z = 14mm: 
% 
% With the mni2tal.m function above on your path, you could 
% type the following at the matlab prompt: 
% 
% 
% mni2tal([10 12 14])
% 
% Which would give the following output (see above): 
% 
% 
% ans =
% 
%     9.9000   12.2692   12.2821
% 
% 
% which is, again, an estimate of the equivalent X, Y and Z 
% coordinates in the Talairach brain. 
% 
% The inverse function, tal2mni.m, gives MNI coordinates for 
% given Talairach coordinates, using the same algorithm. 
% 
% We could of course do a more complex transform to attempt 
% to make a closer match between the two brains. The approach 
% above is only intended to be preliminary. It does have the 
% advantage that it is very simple, and therefore the distortions 
% involved are easy to visualise, and unlikely to have dramatic 
% unexpected effects. 
% 
% Incidentally, if you use the above transform, and you want to 
% cite it, I suggest that you cite this web address. The transform 
% is also mentioned briefly in the following papers: Duncan, J., 
% Seitz, R.J., Kolodny, J., Bor, D., Herzog, H., Ahmed, A., Newell, F.N.,
% Emslie, H. "A neural basis for General Intelligence", Science (21 July
% 2000), 289 (5478), 457-460; Calder, A.J., Lawrence, A.D. and 
% Young,A.W. "Neuropsychology of Fear and Loathing" Nature Reviews 
% Neuroscience (2001), Vol.2 No.5 352-363 
% 
